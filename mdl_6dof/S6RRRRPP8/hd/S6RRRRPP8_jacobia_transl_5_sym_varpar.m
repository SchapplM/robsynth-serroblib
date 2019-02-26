% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function Ja_transl = S6RRRRPP8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:34
% EndTime: 2019-02-26 22:29:34
% DurationCPUTime: 0.20s
% Computational Cost: add. (242->58), mult. (613->99), div. (0->0), fcn. (787->10), ass. (0->39)
t35 = sin(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(2));
t40 = cos(qJ(1));
t47 = cos(pkin(6));
t44 = t40 * t47;
t26 = t35 * t44 + t36 * t39;
t34 = sin(qJ(3));
t38 = cos(qJ(3));
t32 = sin(pkin(6));
t52 = t32 * t40;
t14 = t26 * t38 - t34 * t52;
t25 = t36 * t35 - t39 * t44;
t33 = sin(qJ(4));
t37 = cos(qJ(4));
t1 = t14 * t33 - t25 * t37;
t61 = t14 * t37 + t25 * t33;
t58 = pkin(10) + r_i_i_C(2);
t60 = pkin(3) * t38 + t58 * t34 + pkin(2);
t48 = r_i_i_C(3) + qJ(5);
t59 = pkin(4) + r_i_i_C(1);
t41 = t48 * t33 + t59 * t37 + pkin(3);
t55 = t32 * t36;
t54 = t32 * t38;
t53 = t32 * t39;
t51 = t33 * t38;
t50 = t37 * t38;
t49 = t38 * t39;
t45 = t36 * t47;
t43 = -t26 * t34 - t38 * t52;
t28 = -t35 * t45 + t40 * t39;
t27 = t40 * t35 + t39 * t45;
t24 = t47 * t34 + t35 * t54;
t18 = t28 * t38 + t34 * t55;
t17 = t28 * t34 - t36 * t54;
t11 = t24 * t33 + t37 * t53;
t6 = t18 * t37 + t27 * t33;
t5 = t18 * t33 - t27 * t37;
t2 = [-t36 * pkin(1) - t26 * pkin(2) - t14 * pkin(3) + pkin(8) * t52 - t25 * pkin(9) - t48 * t1 + t58 * t43 - t59 * t61, t28 * pkin(9) + t48 * (-t27 * t51 - t28 * t37) + t59 * (-t27 * t50 + t28 * t33) - t60 * t27, -t41 * t17 + t58 * t18, t48 * t6 - t59 * t5, t5, 0; t40 * pkin(1) + t28 * pkin(2) + t18 * pkin(3) + pkin(8) * t55 + t27 * pkin(9) + t58 * t17 + t48 * t5 + t59 * t6, t26 * pkin(9) + t59 * (-t25 * t50 + t26 * t33) + t48 * (-t25 * t51 - t26 * t37) - t60 * t25, t58 * t14 + t41 * t43, -t59 * t1 + t48 * t61, t1, 0; 0 (t59 * (t33 * t35 + t37 * t49) + t48 * (t33 * t49 - t35 * t37) + pkin(9) * t35 + t60 * t39) * t32, t58 * t24 + t41 * (-t32 * t35 * t34 + t47 * t38) t48 * (t24 * t37 - t33 * t53) - t59 * t11, t11, 0;];
Ja_transl  = t2;
