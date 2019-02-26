% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:51
% EndTime: 2019-02-26 20:23:52
% DurationCPUTime: 0.08s
% Computational Cost: add. (77->28), mult. (154->44), div. (0->0), fcn. (134->7), ass. (0->18)
t52 = pkin(10) + qJ(5);
t50 = sin(t52);
t51 = cos(t52);
t58 = (r_i_i_C(1) * t50 + r_i_i_C(2) * t51) * qJD(5);
t64 = -pkin(1) - pkin(2);
t63 = -r_i_i_C(3) - pkin(7) - qJ(4);
t62 = qJD(5) * t50;
t61 = qJD(5) * t51;
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t56 = sin(qJ(1));
t57 = cos(qJ(1));
t46 = t53 * t57 - t56 * t54;
t47 = t56 * t53 + t54 * t57;
t59 = r_i_i_C(1) * t51 - r_i_i_C(2) * t50 + cos(pkin(10)) * pkin(4) + pkin(3);
t45 = t46 * qJD(1);
t44 = t47 * qJD(1);
t1 = [qJD(2) * t57 - t47 * qJD(4) + t63 * t45 - t46 * t58 - t59 * t44 + (-qJ(2) * t56 + t64 * t57) * qJD(1), qJD(1) * t57, 0, -t44 (-t45 * t51 + t47 * t62) * r_i_i_C(2) + (-t45 * t50 - t47 * t61) * r_i_i_C(1), 0; t56 * qJD(2) + t46 * qJD(4) + t63 * t44 - t47 * t58 + t59 * t45 + (qJ(2) * t57 + t64 * t56) * qJD(1), qJD(1) * t56, 0, t45 (-t44 * t51 - t46 * t62) * r_i_i_C(2) + (-t44 * t50 + t46 * t61) * r_i_i_C(1), 0; 0, 0, 0, 0, t58, 0;];
JaD_transl  = t1;
