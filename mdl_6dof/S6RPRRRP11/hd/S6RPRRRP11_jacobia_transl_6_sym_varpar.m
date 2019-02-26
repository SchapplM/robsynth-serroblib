% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRP11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:31
% EndTime: 2019-02-26 21:13:32
% DurationCPUTime: 0.33s
% Computational Cost: add. (505->66), mult. (1362->110), div. (0->0), fcn. (1790->14), ass. (0->52)
t40 = cos(qJ(1));
t68 = cos(pkin(12));
t70 = cos(pkin(6));
t57 = t70 * t68;
t66 = sin(pkin(12));
t72 = sin(qJ(1));
t50 = -t40 * t57 + t72 * t66;
t33 = sin(pkin(6));
t67 = sin(pkin(7));
t63 = t33 * t67;
t69 = cos(pkin(7));
t82 = t40 * t63 + t50 * t69;
t55 = t70 * t66;
t25 = t40 * t55 + t72 * t68;
t37 = sin(qJ(3));
t73 = cos(qJ(3));
t14 = -t25 * t73 + t82 * t37;
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t64 = t33 * t69;
t42 = -t40 * t64 + t50 * t67;
t81 = t14 * t36 + t42 * t39;
t80 = t14 * t39 - t42 * t36;
t11 = t25 * t37 + t82 * t73;
t77 = pkin(5) + r_i_i_C(1);
t46 = t40 * t66 + t72 * t57;
t75 = t46 * t69 - t72 * t63;
t74 = r_i_i_C(3) + qJ(6) + pkin(11);
t71 = t40 * t33;
t65 = t72 * t33;
t26 = t40 * t68 - t72 * t55;
t15 = t26 * t37 + t75 * t73;
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t16 = t26 * t73 - t75 * t37;
t41 = t46 * t67 + t72 * t64;
t8 = t16 * t39 + t41 * t36;
t1 = t15 * t38 - t8 * t35;
t56 = t70 * t67;
t54 = t69 * t68;
t32 = t38 * pkin(5) + pkin(4);
t53 = t38 * r_i_i_C(1) - t35 * r_i_i_C(2) + t32;
t51 = t38 * r_i_i_C(2) + t77 * t35 + pkin(10);
t45 = -t74 * t36 - t53 * t39 - pkin(3);
t44 = -t68 * t63 + t70 * t69;
t20 = t37 * t56 + (t37 * t54 + t66 * t73) * t33;
t19 = -t73 * t56 + (t37 * t66 - t54 * t73) * t33;
t10 = t20 * t39 + t44 * t36;
t9 = t20 * t36 - t44 * t39;
t7 = t16 * t36 - t41 * t39;
t2 = t15 * t35 + t8 * t38;
t3 = [-t72 * pkin(1) - t25 * pkin(2) + t14 * pkin(3) - t42 * pkin(9) + qJ(2) * t71 - t51 * t11 + t53 * t80 + t74 * t81, t65, t45 * t15 + t51 * t16, -t53 * t7 + t74 * t8, -t2 * r_i_i_C(2) + t77 * t1, t7; qJ(2) * t65 + t40 * pkin(1) + t26 * pkin(2) + t16 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t8 * t32 + t74 * t7 + (t35 * pkin(5) + pkin(10)) * t15 + t41 * pkin(9), -t71, t45 * t11 - t14 * t51, t53 * t81 - t74 * t80 (-t11 * t35 + t38 * t80) * r_i_i_C(2) + t77 * (t11 * t38 + t35 * t80) -t81; 0, t70, t45 * t19 + t51 * t20, t74 * t10 - t53 * t9 (-t10 * t38 - t19 * t35) * r_i_i_C(2) + t77 * (-t10 * t35 + t19 * t38) t9;];
Ja_transl  = t3;
