% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:24
% EndTime: 2019-02-26 22:14:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (354->61), mult. (659->101), div. (0->0), fcn. (845->12), ass. (0->43)
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t44 = cos(qJ(2));
t52 = cos(pkin(6));
t62 = cos(qJ(1));
t47 = t52 * t62;
t26 = t41 * t47 + t42 * t44;
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t38 = sin(pkin(6));
t50 = t38 * t62;
t14 = t26 * t43 - t40 * t50;
t25 = t42 * t41 - t44 * t47;
t36 = pkin(11) + qJ(5);
t34 = sin(t36);
t35 = cos(t36);
t1 = t14 * t34 - t25 * t35;
t66 = t14 * t35 + t25 * t34;
t33 = cos(pkin(11)) * pkin(4) + pkin(3);
t63 = r_i_i_C(2) + pkin(10) + qJ(4);
t65 = t33 * t43 + t63 * t40 + pkin(2);
t64 = pkin(5) + r_i_i_C(1);
t59 = t34 * t43;
t58 = t35 * t43;
t57 = t38 * t42;
t56 = t38 * t43;
t55 = t38 * t44;
t54 = t43 * t44;
t53 = r_i_i_C(3) + qJ(6);
t51 = pkin(4) * sin(pkin(11)) + pkin(9);
t48 = t42 * t52;
t13 = t26 * t40 + t43 * t50;
t45 = -t53 * t34 - t64 * t35 - t33;
t28 = -t41 * t48 + t62 * t44;
t27 = t62 * t41 + t44 * t48;
t24 = t52 * t40 + t41 * t56;
t23 = t38 * t41 * t40 - t52 * t43;
t18 = t28 * t43 + t40 * t57;
t17 = t28 * t40 - t42 * t56;
t11 = t24 * t34 + t35 * t55;
t6 = t18 * t35 + t27 * t34;
t5 = t18 * t34 - t27 * t35;
t2 = [-t42 * pkin(1) - t26 * pkin(2) + pkin(8) * t50 - t53 * t1 - t63 * t13 - t14 * t33 - t51 * t25 - t64 * t66, t53 * (-t27 * t59 - t28 * t35) + t51 * t28 + t64 * (-t27 * t58 + t28 * t34) - t65 * t27, t45 * t17 + t63 * t18, t17, -t64 * t5 + t53 * t6, t5; t62 * pkin(1) + t28 * pkin(2) + pkin(8) * t57 + t63 * t17 + t18 * t33 + t51 * t27 + t53 * t5 + t64 * t6, t64 * (-t25 * t58 + t26 * t34) + t53 * (-t25 * t59 - t26 * t35) + t51 * t26 - t65 * t25, t45 * t13 + t63 * t14, t13, -t64 * t1 + t53 * t66, t1; 0 (t64 * (t34 * t41 + t35 * t54) + t53 * (t34 * t54 - t35 * t41) + t51 * t41 + t65 * t44) * t38, t45 * t23 + t63 * t24, t23, t53 * (t24 * t35 - t34 * t55) - t64 * t11, t11;];
Ja_transl  = t2;
