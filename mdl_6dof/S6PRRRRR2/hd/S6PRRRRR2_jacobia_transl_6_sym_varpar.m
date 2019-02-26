% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:19:14
% EndTime: 2019-02-26 20:19:14
% DurationCPUTime: 0.19s
% Computational Cost: add. (365->50), mult. (507->85), div. (0->0), fcn. (632->14), ass. (0->41)
t71 = r_i_i_C(3) + pkin(11) + pkin(10);
t36 = qJ(5) + qJ(6);
t32 = sin(t36);
t34 = cos(t36);
t45 = cos(qJ(5));
t70 = t45 * pkin(5) + r_i_i_C(1) * t34 - r_i_i_C(2) * t32 + pkin(4);
t37 = qJ(3) + qJ(4);
t33 = sin(t37);
t35 = cos(t37);
t46 = cos(qJ(3));
t69 = t46 * pkin(3) + t71 * t33 + t70 * t35 + pkin(2);
t38 = sin(pkin(12));
t40 = cos(pkin(12));
t47 = cos(qJ(2));
t41 = cos(pkin(6));
t44 = sin(qJ(2));
t57 = t41 * t44;
t26 = t38 * t47 + t40 * t57;
t39 = sin(pkin(6));
t61 = t39 * t40;
t16 = t26 * t35 - t33 * t61;
t56 = t41 * t47;
t25 = t38 * t44 - t40 * t56;
t68 = (-t16 * t32 + t25 * t34) * r_i_i_C(1) + (-t16 * t34 - t25 * t32) * r_i_i_C(2);
t28 = -t38 * t57 + t40 * t47;
t62 = t38 * t39;
t18 = t28 * t35 + t33 * t62;
t27 = t38 * t56 + t40 * t44;
t67 = (-t18 * t32 + t27 * t34) * r_i_i_C(1) + (-t18 * t34 - t27 * t32) * r_i_i_C(2);
t60 = t39 * t44;
t24 = t41 * t33 + t35 * t60;
t58 = t39 * t47;
t63 = (-t24 * t32 - t34 * t58) * r_i_i_C(1) + (-t24 * t34 + t32 * t58) * r_i_i_C(2);
t59 = t39 * t46;
t54 = t71 * t16 + t70 * (-t26 * t33 - t35 * t61);
t53 = t71 * t18 + t70 * (-t28 * t33 + t35 * t62);
t52 = t71 * t24 + t70 * (-t33 * t60 + t41 * t35);
t42 = sin(qJ(5));
t51 = t42 * pkin(5) + t32 * r_i_i_C(1) + t34 * r_i_i_C(2) + pkin(8) + pkin(9);
t43 = sin(qJ(3));
t1 = [0, -t27 * t69 + t51 * t28 (-t28 * t43 + t38 * t59) * pkin(3) + t53, t53 (-t18 * t42 + t27 * t45) * pkin(5) + t67, t67; 0, -t25 * t69 + t51 * t26 (-t26 * t43 - t40 * t59) * pkin(3) + t54, t54 (-t16 * t42 + t25 * t45) * pkin(5) + t68, t68; 1 (t51 * t44 + t69 * t47) * t39 (t41 * t46 - t43 * t60) * pkin(3) + t52, t52 (-t24 * t42 - t45 * t58) * pkin(5) + t63, t63;];
Ja_transl  = t1;
