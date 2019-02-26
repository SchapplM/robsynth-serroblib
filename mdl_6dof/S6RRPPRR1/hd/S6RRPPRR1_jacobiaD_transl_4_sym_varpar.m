% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:07
% EndTime: 2019-02-26 21:28:07
% DurationCPUTime: 0.08s
% Computational Cost: add. (103->20), mult. (160->28), div. (0->0), fcn. (111->6), ass. (0->17)
t145 = qJ(2) + pkin(10);
t143 = sin(t145);
t144 = cos(t145);
t161 = r_i_i_C(3) + qJ(4);
t165 = pkin(3) + r_i_i_C(1);
t154 = t165 * t143 - t161 * t144 + sin(qJ(2)) * pkin(2);
t151 = -t154 * qJD(2) + t143 * qJD(4);
t168 = t151 + (r_i_i_C(2) + qJ(3) + pkin(7)) * qJD(1);
t167 = -t161 * t143 - t165 * t144 - cos(qJ(2)) * pkin(2);
t148 = sin(qJ(1));
t160 = qJD(1) * t148;
t150 = cos(qJ(1));
t159 = qJD(1) * t150;
t158 = qJD(2) * t144;
t153 = qJD(3) + (-pkin(1) + t167) * qJD(1);
t152 = t167 * qJD(2) + qJD(4) * t144;
t1 = [-t168 * t148 + t153 * t150, t152 * t150 + t154 * t160, t159, -t143 * t160 + t150 * t158, 0, 0; t153 * t148 + t168 * t150, t152 * t148 - t154 * t159, t160, t143 * t159 + t148 * t158, 0, 0; 0, t151, 0, qJD(2) * t143, 0, 0;];
JaD_transl  = t1;
