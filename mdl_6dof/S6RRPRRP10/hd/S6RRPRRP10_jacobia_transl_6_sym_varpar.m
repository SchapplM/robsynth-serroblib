% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:03
% EndTime: 2019-02-26 21:51:03
% DurationCPUTime: 0.24s
% Computational Cost: add. (378->63), mult. (631->102), div. (0->0), fcn. (810->12), ass. (0->43)
t40 = sin(qJ(2));
t41 = sin(qJ(1));
t43 = cos(qJ(2));
t44 = cos(qJ(1));
t52 = cos(pkin(6));
t49 = t44 * t52;
t26 = t40 * t49 + t41 * t43;
t35 = pkin(11) + qJ(4);
t33 = sin(t35);
t34 = cos(t35);
t37 = sin(pkin(6));
t56 = t37 * t44;
t14 = t26 * t34 - t33 * t56;
t25 = t41 * t40 - t43 * t49;
t39 = sin(qJ(5));
t42 = cos(qJ(5));
t1 = t14 * t39 - t25 * t42;
t66 = t14 * t42 + t25 * t39;
t32 = cos(pkin(11)) * pkin(3) + pkin(2);
t63 = pkin(10) + r_i_i_C(2);
t65 = pkin(4) * t34 + t63 * t33 + t32;
t53 = r_i_i_C(3) + qJ(6);
t64 = pkin(5) + r_i_i_C(1);
t45 = t53 * t39 + t64 * t42 + pkin(4);
t60 = t34 * t39;
t59 = t34 * t42;
t58 = t37 * t40;
t57 = t37 * t41;
t55 = t39 * t43;
t54 = t42 * t43;
t50 = t41 * t52;
t48 = t37 * (pkin(3) * sin(pkin(11)) + pkin(8));
t47 = -t26 * t33 - t34 * t56;
t38 = -pkin(9) - qJ(3);
t28 = -t40 * t50 + t44 * t43;
t27 = t44 * t40 + t43 * t50;
t22 = t52 * t33 + t34 * t58;
t18 = t28 * t34 + t33 * t57;
t17 = t28 * t33 - t34 * t57;
t11 = t22 * t39 + t37 * t54;
t6 = t18 * t42 + t27 * t39;
t5 = t18 * t39 - t27 * t42;
t2 = [-t41 * pkin(1) - t14 * pkin(4) - t53 * t1 + t25 * t38 - t26 * t32 + t44 * t48 + t63 * t47 - t64 * t66, -t28 * t38 + t53 * (-t27 * t60 - t28 * t42) + t64 * (-t27 * t59 + t28 * t39) - t65 * t27, t27, -t45 * t17 + t63 * t18, -t64 * t5 + t53 * t6, t5; t44 * pkin(1) + t18 * pkin(4) + t63 * t17 - t27 * t38 + t28 * t32 + t41 * t48 + t53 * t5 + t64 * t6, -t26 * t38 + t64 * (-t25 * t59 + t26 * t39) + t53 * (-t25 * t60 - t26 * t42) - t65 * t25, t25, t63 * t14 + t45 * t47, -t64 * t1 + t53 * t66, t1; 0 (t64 * (t34 * t54 + t39 * t40) + t53 * (t34 * t55 - t40 * t42) - t38 * t40 + t65 * t43) * t37, -t37 * t43, t63 * t22 + t45 * (-t33 * t58 + t34 * t52) t53 * (t22 * t42 - t37 * t55) - t64 * t11, t11;];
Ja_transl  = t2;
