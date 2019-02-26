% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PPRRRP2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRRRP2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:04
% EndTime: 2019-02-26 19:42:05
% DurationCPUTime: 0.23s
% Computational Cost: add. (475->58), mult. (1331->105), div. (0->0), fcn. (1760->14), ass. (0->52)
t39 = cos(pkin(11));
t61 = sin(pkin(12));
t62 = sin(pkin(11));
t52 = t62 * t61;
t64 = cos(pkin(12));
t59 = t39 * t64;
t66 = cos(pkin(6));
t47 = -t66 * t59 + t52;
t38 = sin(pkin(6));
t63 = sin(pkin(7));
t60 = t38 * t63;
t65 = cos(pkin(7));
t75 = t39 * t60 + t47 * t65;
t53 = t62 * t64;
t58 = t39 * t61;
t48 = t66 * t53 + t58;
t57 = t62 * t38;
t74 = t48 * t65 - t63 * t57;
t73 = pkin(5) + r_i_i_C(1);
t72 = pkin(10) + r_i_i_C(2);
t71 = cos(qJ(3));
t70 = t39 * t38;
t40 = sin(qJ(5));
t44 = cos(qJ(4));
t69 = t40 * t44;
t43 = cos(qJ(5));
t68 = t43 * t44;
t67 = r_i_i_C(3) + qJ(6);
t55 = t65 * t64;
t54 = t63 * t66;
t41 = sin(qJ(4));
t50 = -pkin(4) * t44 - t72 * t41 - pkin(3);
t49 = t67 * t40 + t73 * t43 + pkin(4);
t42 = sin(qJ(3));
t34 = -t66 * t52 + t59;
t33 = t66 * t58 + t53;
t32 = -t64 * t60 + t66 * t65;
t29 = t48 * t63 + t65 * t57;
t28 = t47 * t63 - t65 * t70;
t27 = t42 * t54 + (t42 * t55 + t71 * t61) * t38;
t26 = -t71 * t54 + (t42 * t61 - t55 * t71) * t38;
t24 = t27 * t44 + t32 * t41;
t22 = t34 * t71 - t74 * t42;
t21 = t34 * t42 + t74 * t71;
t20 = t33 * t71 - t75 * t42;
t19 = t33 * t42 + t75 * t71;
t13 = t24 * t40 - t26 * t43;
t12 = t22 * t44 + t29 * t41;
t10 = t20 * t44 + t28 * t41;
t3 = t12 * t40 - t21 * t43;
t1 = t10 * t40 - t19 * t43;
t2 = [0, t57, t22 * pkin(9) + t73 * (-t21 * t68 + t22 * t40) + t67 * (-t21 * t69 - t22 * t43) + t50 * t21, t72 * t12 + t49 * (-t22 * t41 + t29 * t44) t67 * (t12 * t43 + t21 * t40) - t73 * t3, t3; 0, -t70, t20 * pkin(9) + t73 * (-t19 * t68 + t20 * t40) + t67 * (-t19 * t69 - t20 * t43) + t50 * t19, t72 * t10 + t49 * (-t20 * t41 + t28 * t44) t67 * (t10 * t43 + t19 * t40) - t73 * t1, t1; 1, t66, t27 * pkin(9) + t73 * (-t26 * t68 + t27 * t40) + t67 * (-t26 * t69 - t27 * t43) + t50 * t26, t72 * t24 + t49 * (-t27 * t41 + t32 * t44) t67 * (t24 * t43 + t26 * t40) - t73 * t13, t13;];
Ja_transl  = t2;
