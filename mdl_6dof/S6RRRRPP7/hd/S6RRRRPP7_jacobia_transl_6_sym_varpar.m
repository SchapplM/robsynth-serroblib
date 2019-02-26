% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP7
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:54
% EndTime: 2019-02-26 22:28:55
% DurationCPUTime: 0.26s
% Computational Cost: add. (367->67), mult. (693->110), div. (0->0), fcn. (886->12), ass. (0->45)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t45 = cos(qJ(2));
t53 = cos(pkin(6));
t63 = cos(qJ(1));
t48 = t53 * t63;
t26 = t41 * t48 + t42 * t45;
t40 = sin(qJ(3));
t44 = cos(qJ(3));
t37 = sin(pkin(6));
t51 = t37 * t63;
t14 = t26 * t44 - t40 * t51;
t25 = t42 * t41 - t45 * t48;
t36 = qJ(4) + pkin(11);
t34 = sin(t36);
t35 = cos(t36);
t1 = t14 * t34 - t25 * t35;
t67 = t14 * t35 + t25 * t34;
t43 = cos(qJ(4));
t33 = t43 * pkin(4) + pkin(3);
t64 = r_i_i_C(2) + qJ(5) + pkin(10);
t66 = t33 * t44 + t64 * t40 + pkin(2);
t65 = pkin(5) + r_i_i_C(1);
t60 = t34 * t44;
t59 = t35 * t44;
t58 = t37 * t42;
t57 = t37 * t44;
t56 = t37 * t45;
t55 = t44 * t45;
t54 = r_i_i_C(3) + qJ(6);
t39 = sin(qJ(4));
t52 = pkin(4) * t39 + pkin(9);
t49 = t42 * t53;
t13 = t26 * t40 + t44 * t51;
t46 = -t54 * t34 - t65 * t35 - t33;
t28 = -t41 * t49 + t63 * t45;
t27 = t63 * t41 + t45 * t49;
t24 = t53 * t40 + t41 * t57;
t23 = t37 * t41 * t40 - t53 * t44;
t18 = t28 * t44 + t40 * t58;
t17 = t28 * t40 - t42 * t57;
t11 = t24 * t34 + t35 * t56;
t6 = t18 * t35 + t27 * t34;
t5 = t18 * t34 - t27 * t35;
t2 = [-t42 * pkin(1) - t26 * pkin(2) + pkin(8) * t51 - t54 * t1 - t64 * t13 - t14 * t33 - t52 * t25 - t65 * t67, t54 * (-t27 * t60 - t28 * t35) + t52 * t28 + t65 * (-t27 * t59 + t28 * t34) - t66 * t27, t46 * t17 + t64 * t18, t54 * t6 - t65 * t5 + (-t18 * t39 + t27 * t43) * pkin(4), t17, t5; t63 * pkin(1) + t28 * pkin(2) + pkin(8) * t58 + t64 * t17 + t18 * t33 + t52 * t27 + t54 * t5 + t65 * t6, t65 * (-t25 * t59 + t26 * t34) + t54 * (-t25 * t60 - t26 * t35) + t52 * t26 - t66 * t25, t46 * t13 + t64 * t14, t54 * t67 - t65 * t1 + (-t14 * t39 + t25 * t43) * pkin(4), t13, t1; 0 (t65 * (t34 * t41 + t35 * t55) + t54 * (t34 * t55 - t35 * t41) + t52 * t41 + t66 * t45) * t37, t46 * t23 + t64 * t24, t54 * (t24 * t35 - t34 * t56) - t65 * t11 + (-t24 * t39 - t43 * t56) * pkin(4), t23, t11;];
Ja_transl  = t2;
