% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:36
% EndTime: 2019-02-26 20:28:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (53->25), mult. (150->38), div. (0->0), fcn. (104->4), ass. (0->16)
t139 = sin(qJ(4));
t141 = cos(qJ(4));
t151 = r_i_i_C(3) + qJ(5);
t152 = pkin(4) - r_i_i_C(2);
t143 = -(t151 * t139 + t152 * t141) * qJD(4) + t141 * qJD(5);
t154 = -qJD(3) + t143;
t140 = sin(qJ(1));
t150 = qJD(1) * t140;
t142 = cos(qJ(1));
t138 = qJD(1) * t142;
t149 = qJD(4) * t139;
t147 = pkin(7) + r_i_i_C(1) - qJ(2);
t146 = qJD(4) * t151;
t145 = -t152 * qJD(4) + qJD(5);
t144 = -t152 * t139 + t151 * t141 - pkin(1) - qJ(3);
t1 = [t142 * qJD(2) + t154 * t140 + (t147 * t140 + t144 * t142) * qJD(1), t138, -t150 (t142 * t146 - t152 * t150) * t141 + (t145 * t142 - t151 * t150) * t139, t141 * t150 + t142 * t149, 0; t140 * qJD(2) - t154 * t142 + (t144 * t140 - t147 * t142) * qJD(1), t150, t138 (t152 * t138 + t140 * t146) * t141 + (t151 * t138 + t145 * t140) * t139, -t141 * t138 + t140 * t149, 0; 0, 0, 0, t143, qJD(4) * t141, 0;];
JaD_transl  = t1;
