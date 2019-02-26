% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:20
% EndTime: 2019-02-26 20:40:20
% DurationCPUTime: 0.08s
% Computational Cost: add. (92->20), mult. (138->30), div. (0->0), fcn. (94->6), ass. (0->17)
t144 = sin(qJ(3));
t145 = cos(qJ(3));
t155 = r_i_i_C(3) + qJ(4);
t157 = pkin(3) + r_i_i_C(1);
t150 = t157 * t144 - t155 * t145;
t159 = qJD(1) * t150;
t158 = t150 * qJD(3) - t144 * qJD(4);
t156 = pkin(7) + r_i_i_C(2);
t154 = qJD(1) * t144;
t153 = qJD(3) * t145;
t151 = -t155 * t144 - t157 * t145;
t148 = -pkin(2) + t151;
t147 = t151 * qJD(3) + qJD(4) * t145;
t143 = qJ(1) + pkin(9);
t142 = cos(t143);
t141 = sin(t143);
t1 = [t158 * t141 + (-cos(qJ(1)) * pkin(1) - t156 * t141 + t148 * t142) * qJD(1), 0, t141 * t159 + t147 * t142, -t141 * t154 + t142 * t153, 0, 0; -t158 * t142 + (-sin(qJ(1)) * pkin(1) + t156 * t142 + t148 * t141) * qJD(1), 0, t147 * t141 - t142 * t159, t141 * t153 + t142 * t154, 0, 0; 0, 0, -t158, qJD(3) * t144, 0, 0;];
JaD_transl  = t1;
