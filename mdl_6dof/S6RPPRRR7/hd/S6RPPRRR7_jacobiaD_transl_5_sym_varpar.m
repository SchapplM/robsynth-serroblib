% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:38:12
% EndTime: 2019-02-26 20:38:12
% DurationCPUTime: 0.09s
% Computational Cost: add. (129->32), mult. (132->44), div. (0->0), fcn. (87->7), ass. (0->28)
t48 = pkin(10) + qJ(4);
t44 = cos(t48);
t67 = pkin(4) * t44;
t45 = qJ(5) + t48;
t42 = cos(t45);
t66 = r_i_i_C(1) * t42;
t41 = sin(t45);
t65 = r_i_i_C(2) * t41;
t49 = qJD(4) + qJD(5);
t64 = t41 * t49;
t63 = t42 * t49;
t50 = sin(qJ(1));
t62 = qJD(1) * t50;
t51 = cos(qJ(1));
t46 = qJD(1) * t51;
t43 = sin(t48);
t61 = qJD(4) * t43;
t60 = -pkin(1) - r_i_i_C(3) - pkin(8) - pkin(7) - qJ(3);
t59 = r_i_i_C(1) * t64;
t58 = qJD(4) * t67;
t57 = qJD(1) * t66;
t56 = -r_i_i_C(1) * t63 + r_i_i_C(2) * t64;
t55 = -r_i_i_C(1) * t41 - r_i_i_C(2) * t42;
t54 = qJ(2) + pkin(4) * t43 + sin(pkin(10)) * pkin(3) - t55;
t53 = t50 * t57 - t62 * t65 + (t63 * r_i_i_C(2) + t59) * t51;
t52 = t58 + qJD(2) + (-t65 + t66) * t49;
t39 = t51 * t57;
t1 = [-t50 * qJD(3) + t52 * t51 + (-t54 * t50 + t60 * t51) * qJD(1), t46, -t62, t39 + (-t65 + t67) * t46 + (-pkin(4) * t61 + t55 * t49) * t50, -t50 * t59 + t39 + (-t41 * t46 - t50 * t63) * r_i_i_C(2), 0; t51 * qJD(3) + t52 * t50 + (t60 * t50 + t54 * t51) * qJD(1), t62, t46 (t44 * t62 + t51 * t61) * pkin(4) + t53, t53, 0; 0, 0, 0, t56 - t58, t56, 0;];
JaD_transl  = t1;
