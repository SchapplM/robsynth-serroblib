% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:45
% EndTime: 2019-02-26 20:50:45
% DurationCPUTime: 0.08s
% Computational Cost: add. (92->20), mult. (138->30), div. (0->0), fcn. (94->6), ass. (0->17)
t147 = sin(qJ(3));
t148 = cos(qJ(3));
t158 = r_i_i_C(3) + qJ(4);
t160 = pkin(3) - r_i_i_C(2);
t153 = t160 * t147 - t158 * t148;
t162 = qJD(1) * t153;
t161 = t153 * qJD(3) - t147 * qJD(4);
t159 = pkin(7) + r_i_i_C(1);
t157 = qJD(1) * t147;
t156 = qJD(3) * t148;
t154 = -t158 * t147 - t160 * t148;
t151 = -pkin(2) + t154;
t150 = t154 * qJD(3) + qJD(4) * t148;
t146 = qJ(1) + pkin(10);
t145 = cos(t146);
t144 = sin(t146);
t1 = [t161 * t144 + (-cos(qJ(1)) * pkin(1) - t159 * t144 + t151 * t145) * qJD(1), 0, t144 * t162 + t150 * t145, -t144 * t157 + t145 * t156, 0, 0; -t161 * t145 + (-sin(qJ(1)) * pkin(1) + t159 * t145 + t151 * t144) * qJD(1), 0, t150 * t144 - t145 * t162, t144 * t156 + t145 * t157, 0, 0; 0, 0, -t161, qJD(3) * t147, 0, 0;];
JaD_transl  = t1;
