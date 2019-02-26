% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_jacobia_transl_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:07
% EndTime: 2019-02-26 21:21:09
% DurationCPUTime: 0.53s
% Computational Cost: add. (654->92), mult. (1858->160), div. (0->0), fcn. (2457->16), ass. (0->66)
t47 = cos(pkin(6));
t44 = cos(pkin(14));
t55 = cos(qJ(1));
t66 = t55 * t44;
t40 = sin(pkin(14));
t51 = sin(qJ(1));
t71 = t51 * t40;
t34 = -t47 * t66 + t71;
t68 = t55 * t40;
t69 = t51 * t44;
t35 = t47 * t68 + t69;
t46 = cos(pkin(7));
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t42 = sin(pkin(7));
t43 = sin(pkin(6));
t67 = t55 * t43;
t64 = t42 * t67;
t24 = (t34 * t46 + t64) * t54 + t35 * t50;
t31 = -t34 * t42 + t46 * t67;
t41 = sin(pkin(8));
t45 = cos(pkin(8));
t15 = t24 * t41 - t31 * t45;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t72 = t46 * t50;
t25 = t34 * t72 - t35 * t54 + t50 * t64;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t62 = t24 * t45 + t31 * t41;
t6 = t25 * t53 + t62 * t49;
t89 = t15 * t52 + t6 * t48;
t88 = -t15 * t48 + t6 * t52;
t85 = t25 * t49 - t62 * t53;
t80 = r_i_i_C(3) + pkin(12);
t79 = t41 * pkin(11);
t77 = t41 * t48;
t76 = t41 * t52;
t75 = t42 * t47;
t74 = t45 * t49;
t73 = t45 * t53;
t70 = t51 * t43;
t65 = t43 * qJ(2);
t29 = t54 * t75 + (t44 * t46 * t54 - t40 * t50) * t43;
t33 = -t43 * t44 * t42 + t47 * t46;
t61 = t29 * t45 + t33 * t41;
t60 = t52 * r_i_i_C(1) - t48 * r_i_i_C(2) + pkin(4);
t36 = -t47 * t69 - t68;
t58 = -t36 * t42 + t46 * t70;
t57 = t36 * t46 + t42 * t70;
t56 = t58 * t41;
t37 = -t47 * t71 + t66;
t26 = -t37 * t50 + t57 * t54;
t17 = -t26 * t41 + t58 * t45;
t30 = t50 * t75 + (t40 * t54 + t44 * t72) * t43;
t27 = t37 * t54 + t57 * t50;
t21 = -t29 * t41 + t33 * t45;
t19 = t29 * t53 - t30 * t74;
t14 = t30 * t53 + t61 * t49;
t12 = t26 * t53 - t27 * t74;
t10 = -t24 * t53 + t25 * t74;
t8 = t27 * t53 + (t26 * t45 + t56) * t49;
t7 = -t26 * t73 + t27 * t49 - t53 * t56;
t2 = t17 * t48 + t8 * t52;
t1 = t17 * t52 - t8 * t48;
t3 = [-t51 * pkin(1) - t35 * pkin(2) + t25 * pkin(3) + t6 * pkin(4) + t31 * pkin(10) - t15 * pkin(11) + t88 * r_i_i_C(1) - t89 * r_i_i_C(2) + t55 * t65 + t80 * t85, t70 (t12 * t52 + t27 * t77) * r_i_i_C(1) + (-t12 * t48 + t27 * t76) * r_i_i_C(2) + t12 * pkin(4) + t26 * pkin(3) + t27 * t79 + t80 * (t26 * t49 + t27 * t73) -t60 * t7 + t80 * t8, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t55 * pkin(1) + t37 * pkin(2) + t27 * pkin(3) + t8 * pkin(4) + t58 * pkin(10) + t17 * pkin(11) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t51 * t65 + t80 * t7, -t67 (t10 * t52 - t25 * t77) * r_i_i_C(1) + (-t10 * t48 - t25 * t76) * r_i_i_C(2) + t10 * pkin(4) - t24 * pkin(3) - t25 * t79 + t80 * (-t24 * t49 - t25 * t73) -t6 * t80 + t60 * t85, t89 * r_i_i_C(1) + t88 * r_i_i_C(2), 0; 0, t47 (t19 * t52 + t30 * t77) * r_i_i_C(1) + (-t19 * t48 + t30 * t76) * r_i_i_C(2) + t19 * pkin(4) + t29 * pkin(3) + t30 * t79 + t80 * (t29 * t49 + t30 * t73) t80 * t14 + t60 * (-t30 * t49 + t61 * t53) (-t14 * t48 + t21 * t52) * r_i_i_C(1) + (-t14 * t52 - t21 * t48) * r_i_i_C(2), 0;];
Ja_transl  = t3;
