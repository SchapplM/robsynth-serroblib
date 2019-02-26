% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:41
% EndTime: 2019-02-26 20:14:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (284->61), mult. (791->116), div. (0->0), fcn. (1027->12), ass. (0->51)
t66 = pkin(4) - r_i_i_C(2);
t65 = pkin(10) + r_i_i_C(1);
t34 = sin(pkin(7));
t64 = pkin(9) * t34;
t35 = sin(pkin(6));
t63 = t34 * t35;
t38 = cos(pkin(6));
t62 = t34 * t38;
t39 = sin(qJ(4));
t61 = t34 * t39;
t41 = sin(qJ(2));
t60 = t34 * t41;
t42 = cos(qJ(4));
t59 = t34 * t42;
t37 = cos(pkin(7));
t58 = t35 * t37;
t40 = sin(qJ(3));
t57 = t37 * t40;
t43 = cos(qJ(3));
t56 = t37 * t43;
t55 = t38 * t41;
t44 = cos(qJ(2));
t54 = t38 * t44;
t53 = t40 * t41;
t52 = t40 * t44;
t51 = t41 * t43;
t50 = t43 * t44;
t49 = r_i_i_C(3) + qJ(5);
t48 = t35 * t60;
t33 = sin(pkin(12));
t36 = cos(pkin(12));
t28 = -t33 * t41 + t36 * t54;
t47 = t28 * t37 - t36 * t63;
t30 = -t33 * t54 - t36 * t41;
t46 = t30 * t37 + t33 * t63;
t45 = t49 * t39 + t66 * t42 + pkin(3);
t31 = -t33 * t55 + t36 * t44;
t29 = t33 * t44 + t36 * t55;
t27 = t38 * t37 - t44 * t63;
t26 = (-t37 * t53 + t50) * t35;
t24 = -t30 * t34 + t33 * t58;
t23 = -t28 * t34 - t36 * t58;
t22 = t40 * t62 + (t37 * t52 + t51) * t35;
t18 = t30 * t43 - t31 * t57;
t16 = t28 * t43 - t29 * t57;
t13 = t22 * t39 - t27 * t42;
t12 = t31 * t43 + t46 * t40;
t10 = t29 * t43 + t47 * t40;
t3 = t12 * t39 - t24 * t42;
t1 = t10 * t39 - t23 * t42;
t2 = [0, t31 * t64 + t30 * pkin(2) + t18 * pkin(3) + t66 * (t18 * t42 + t31 * t61) + t49 * (t18 * t39 - t31 * t59) + t65 * (t30 * t40 + t31 * t56) t65 * t12 + t45 * (-t31 * t40 + t46 * t43) t49 * (t12 * t42 + t24 * t39) - t66 * t3, t3, 0; 0, t29 * t64 + t28 * pkin(2) + t16 * pkin(3) + t66 * (t16 * t42 + t29 * t61) + t49 * (t16 * t39 - t29 * t59) + t65 * (t28 * t40 + t29 * t56) t65 * t10 + t45 * (-t29 * t40 + t47 * t43) t49 * (t10 * t42 + t23 * t39) - t66 * t1, t1, 0; 1, t26 * pkin(3) + t66 * (t26 * t42 + t39 * t48) + t49 * (t26 * t39 - t42 * t48) + (pkin(2) * t44 + pkin(9) * t60 + t65 * (t37 * t51 + t52)) * t35, t65 * t22 + t45 * (t43 * t62 + (t37 * t50 - t53) * t35) t49 * (t22 * t42 + t27 * t39) - t66 * t13, t13, 0;];
Ja_transl  = t2;
