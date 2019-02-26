% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:33
% EndTime: 2019-02-26 22:29:34
% DurationCPUTime: 0.21s
% Computational Cost: add. (307->59), mult. (775->98), div. (0->0), fcn. (998->10), ass. (0->40)
t36 = sin(qJ(2));
t39 = cos(qJ(2));
t49 = cos(pkin(6));
t59 = cos(qJ(1));
t43 = t49 * t59;
t58 = sin(qJ(1));
t26 = t36 * t43 + t58 * t39;
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t33 = sin(pkin(6));
t45 = t33 * t59;
t14 = t26 * t38 - t35 * t45;
t25 = t58 * t36 - t39 * t43;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t1 = t14 * t34 - t25 * t37;
t61 = t14 * t37 + t25 * t34;
t47 = pkin(10) - r_i_i_C(3) - qJ(6);
t60 = pkin(3) * t38 + t47 * t35 + pkin(2);
t48 = pkin(4) + pkin(5) + r_i_i_C(1);
t50 = r_i_i_C(2) + qJ(5);
t40 = t50 * t34 + t48 * t37 + pkin(3);
t55 = t33 * t36;
t54 = t33 * t39;
t53 = t34 * t38;
t52 = t37 * t38;
t51 = t38 * t39;
t44 = t33 * t58;
t13 = -t26 * t35 - t38 * t45;
t42 = t49 * t58;
t28 = -t36 * t42 + t59 * t39;
t27 = t59 * t36 + t39 * t42;
t24 = t49 * t35 + t38 * t55;
t23 = -t35 * t55 + t49 * t38;
t18 = t28 * t38 + t35 * t44;
t17 = t28 * t35 - t38 * t44;
t11 = t24 * t34 + t37 * t54;
t6 = t18 * t37 + t27 * t34;
t5 = t18 * t34 - t27 * t37;
t2 = [-t58 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t45 - t25 * pkin(9) - t50 * t1 + t47 * t13 - t48 * t61, t28 * pkin(9) + t50 * (-t27 * t53 - t28 * t37) + t48 * (-t27 * t52 + t28 * t34) - t60 * t27, -t40 * t17 + t47 * t18, -t48 * t5 + t50 * t6, t5, -t17; t59 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t44 + t27 * pkin(9) + t47 * t17 + t48 * t6 + t50 * t5, t26 * pkin(9) + t50 * (-t25 * t53 - t26 * t37) + t48 * (-t25 * t52 + t26 * t34) - t60 * t25, t40 * t13 + t47 * t14, -t48 * t1 + t50 * t61, t1, t13; 0 (t50 * (t34 * t51 - t36 * t37) + t48 * (t34 * t36 + t37 * t51) + pkin(9) * t36 + t60 * t39) * t33, t40 * t23 + t47 * t24, t50 * (t24 * t37 - t34 * t54) - t48 * t11, t11, t23;];
Ja_transl  = t2;
