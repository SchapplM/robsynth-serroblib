% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR7
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
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:06
% EndTime: 2019-02-26 22:34:06
% DurationCPUTime: 0.22s
% Computational Cost: add. (469->66), mult. (550->99), div. (0->0), fcn. (685->14), ass. (0->46)
t64 = pkin(11) + r_i_i_C(3);
t40 = sin(qJ(6));
t43 = cos(qJ(6));
t68 = t43 * r_i_i_C(1) - t40 * r_i_i_C(2) + pkin(5);
t41 = sin(qJ(2));
t42 = sin(qJ(1));
t44 = cos(qJ(2));
t45 = cos(qJ(1));
t55 = cos(pkin(6));
t52 = t45 * t55;
t21 = t41 * t52 + t42 * t44;
t38 = qJ(3) + qJ(4);
t33 = pkin(12) + t38;
t30 = sin(t33);
t31 = cos(t33);
t39 = sin(pkin(6));
t56 = t39 * t45;
t10 = t21 * t31 - t30 * t56;
t20 = t42 * t41 - t44 * t52;
t67 = t10 * t40 - t20 * t43;
t66 = -t10 * t43 - t20 * t40;
t35 = cos(t38);
t28 = pkin(4) * t35 + cos(qJ(3)) * pkin(3);
t26 = pkin(2) + t28;
t65 = t64 * t30 + t68 * t31 + t26;
t59 = t39 * t41;
t58 = t39 * t42;
t57 = t39 * t44;
t34 = sin(t38);
t27 = pkin(4) * t34 + sin(qJ(3)) * pkin(3);
t54 = t39 * (pkin(8) + t27);
t53 = t42 * t55;
t37 = -qJ(5) - pkin(10) - pkin(9);
t50 = t40 * r_i_i_C(1) + t43 * r_i_i_C(2) - t37;
t9 = -t21 * t30 - t31 * t56;
t49 = t64 * t10 + t68 * t9;
t23 = -t41 * t53 + t45 * t44;
t13 = t23 * t30 - t31 * t58;
t14 = t23 * t31 + t30 * t58;
t48 = -t68 * t13 + t64 * t14;
t19 = t55 * t30 + t31 * t59;
t47 = t64 * t19 + t68 * (-t30 * t59 + t55 * t31);
t22 = t45 * t41 + t44 * t53;
t2 = t14 * t43 + t22 * t40;
t1 = -t14 * t40 + t22 * t43;
t3 = [-t42 * pkin(1) - t10 * pkin(5) + t66 * r_i_i_C(1) + t67 * r_i_i_C(2) + t20 * t37 - t21 * t26 + t45 * t54 + t64 * t9, -t22 * t65 + t50 * t23, -t23 * t27 + t28 * t58 + t48 (-t23 * t34 + t35 * t58) * pkin(4) + t48, t22, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t45 * pkin(1) + t14 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t64 * t13 - t22 * t37 + t23 * t26 + t42 * t54, -t20 * t65 + t50 * t21, -t21 * t27 - t28 * t56 + t49 (-t21 * t34 - t35 * t56) * pkin(4) + t49, t20, -t67 * r_i_i_C(1) + t66 * r_i_i_C(2); 0 (t50 * t41 + t65 * t44) * t39, -t27 * t59 + t55 * t28 + t47 (-t34 * t59 + t55 * t35) * pkin(4) + t47, -t57 (-t19 * t40 - t43 * t57) * r_i_i_C(1) + (-t19 * t43 + t40 * t57) * r_i_i_C(2);];
Ja_transl  = t3;
