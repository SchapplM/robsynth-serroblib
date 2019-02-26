% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RPRRRP11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:13:41
% EndTime: 2019-02-26 21:13:41
% DurationCPUTime: 0.27s
% Computational Cost: add. (414->64), mult. (1142->110), div. (0->0), fcn. (1503->14), ass. (0->51)
t39 = cos(qJ(1));
t56 = sin(pkin(12));
t60 = cos(pkin(6));
t48 = t60 * t56;
t58 = cos(pkin(12));
t62 = sin(qJ(1));
t26 = t39 * t48 + t58 * t62;
t36 = sin(qJ(3));
t63 = cos(qJ(3));
t50 = t60 * t58;
t25 = -t39 * t50 + t56 * t62;
t33 = sin(pkin(6));
t57 = sin(pkin(7));
t52 = t33 * t57;
t59 = cos(pkin(7));
t67 = t25 * t59 + t39 * t52;
t11 = t26 * t36 + t63 * t67;
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t12 = t26 * t63 - t36 * t67;
t53 = t33 * t59;
t21 = t25 * t57 - t39 * t53;
t35 = sin(qJ(4));
t38 = cos(qJ(4));
t4 = t12 * t38 + t21 * t35;
t71 = -t11 * t37 + t34 * t4;
t70 = -t11 * t34 - t37 * t4;
t66 = -t12 * t35 + t21 * t38;
t42 = t39 * t56 + t50 * t62;
t65 = t42 * t59 - t52 * t62;
t64 = r_i_i_C(3) + pkin(11);
t61 = t39 * t33;
t55 = t62 * t33;
t49 = t60 * t57;
t47 = t59 * t58;
t46 = r_i_i_C(1) * t37 - r_i_i_C(2) * t34 + pkin(4);
t45 = r_i_i_C(1) * t34 + r_i_i_C(2) * t37 + pkin(10);
t43 = -t35 * t64 - t38 * t46 - pkin(3);
t40 = t42 * t57 + t53 * t62;
t27 = t39 * t58 - t48 * t62;
t24 = -t52 * t58 + t59 * t60;
t19 = t36 * t49 + (t36 * t47 + t56 * t63) * t33;
t18 = -t63 * t49 + (t36 * t56 - t47 * t63) * t33;
t16 = t27 * t63 - t36 * t65;
t15 = t27 * t36 + t63 * t65;
t10 = t19 * t38 + t24 * t35;
t8 = t16 * t38 + t35 * t40;
t7 = t16 * t35 - t38 * t40;
t2 = t15 * t34 + t37 * t8;
t1 = t15 * t37 - t34 * t8;
t3 = [-t62 * pkin(1) - t26 * pkin(2) - t12 * pkin(3) - pkin(4) * t4 - t21 * pkin(9) - t11 * pkin(10) + r_i_i_C(1) * t70 + r_i_i_C(2) * t71 + qJ(2) * t61 + t64 * t66, t55, t15 * t43 + t16 * t45, -t46 * t7 + t64 * t8, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t39 * pkin(1) + t27 * pkin(2) + t16 * pkin(3) + t8 * pkin(4) + pkin(9) * t40 + t15 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t55 + t64 * t7, -t61, t11 * t43 + t12 * t45, t4 * t64 + t46 * t66, -r_i_i_C(1) * t71 + r_i_i_C(2) * t70, 0; 0, t60, t18 * t43 + t19 * t45, t64 * t10 + t46 * (-t19 * t35 + t24 * t38) (-t10 * t34 + t18 * t37) * r_i_i_C(1) + (-t10 * t37 - t18 * t34) * r_i_i_C(2), 0;];
Ja_transl  = t3;
