% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR9_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR9_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR9_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR9_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:35:12
% EndTime: 2019-02-26 22:35:13
% DurationCPUTime: 0.21s
% Computational Cost: add. (400->57), mult. (558->91), div. (0->0), fcn. (698->14), ass. (0->44)
t64 = r_i_i_C(3) + pkin(11) + qJ(5);
t31 = cos(pkin(12)) * pkin(5) + pkin(4);
t37 = pkin(12) + qJ(6);
t33 = sin(t37);
t34 = cos(t37);
t54 = t34 * r_i_i_C(1) - t33 * r_i_i_C(2) + t31;
t45 = cos(qJ(3));
t32 = t45 * pkin(3) + pkin(2);
t38 = qJ(3) + qJ(4);
t35 = sin(t38);
t36 = cos(t38);
t67 = t64 * t35 + t54 * t36 + t32;
t40 = sin(pkin(6));
t43 = sin(qJ(2));
t63 = t40 * t43;
t44 = sin(qJ(1));
t62 = t40 * t44;
t46 = cos(qJ(2));
t61 = t40 * t46;
t47 = cos(qJ(1));
t60 = t40 * t47;
t59 = cos(pkin(6));
t58 = sin(pkin(12)) * pkin(5) + pkin(10) + pkin(9);
t56 = t47 * t59;
t24 = t43 * t56 + t44 * t46;
t12 = t24 * t36 - t35 * t60;
t57 = t44 * t59;
t42 = sin(qJ(3));
t55 = t40 * (pkin(3) * t42 + pkin(8));
t11 = t24 * t35 + t36 * t60;
t53 = -t54 * t11 + t64 * t12;
t26 = -t43 * t57 + t47 * t46;
t15 = t26 * t35 - t36 * t62;
t16 = t26 * t36 + t35 * t62;
t52 = -t54 * t15 + t64 * t16;
t21 = t35 * t63 - t59 * t36;
t22 = t59 * t35 + t36 * t63;
t51 = -t54 * t21 + t64 * t22;
t50 = t33 * r_i_i_C(1) + t34 * r_i_i_C(2) + t58;
t25 = t47 * t43 + t46 * t57;
t23 = t44 * t43 - t46 * t56;
t2 = t16 * t34 + t25 * t33;
t1 = -t16 * t33 + t25 * t34;
t3 = [-t44 * pkin(1) - t64 * t11 - t54 * t12 - t50 * t23 - t24 * t32 + t47 * t55, -t25 * t67 + t50 * t26 (-t26 * t42 + t45 * t62) * pkin(3) + t52, t52, t15, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t47 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t64 * t15 + t16 * t31 + t58 * t25 + t26 * t32 + t44 * t55, -t23 * t67 + t50 * t24 (-t24 * t42 - t45 * t60) * pkin(3) + t53, t53, t11 (-t12 * t33 + t23 * t34) * r_i_i_C(1) + (-t12 * t34 - t23 * t33) * r_i_i_C(2); 0 (t50 * t43 + t67 * t46) * t40 (-t42 * t63 + t59 * t45) * pkin(3) + t51, t51, t21 (-t22 * t33 - t34 * t61) * r_i_i_C(1) + (-t22 * t34 + t33 * t61) * r_i_i_C(2);];
Ja_transl  = t3;
