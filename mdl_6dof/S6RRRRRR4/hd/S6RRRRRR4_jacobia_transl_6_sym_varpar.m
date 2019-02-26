% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:47
% EndTime: 2019-02-26 22:48:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (345->44), mult. (217->57), div. (0->0), fcn. (237->12), ass. (0->35)
t28 = cos(qJ(2));
t26 = sin(qJ(2));
t39 = r_i_i_C(3) + pkin(11) + pkin(10) + pkin(9) + pkin(8);
t33 = t39 * t26;
t25 = qJ(3) + qJ(4);
t22 = qJ(5) + t25;
t19 = cos(t22);
t13 = pkin(5) * t19 + pkin(4) * cos(t25);
t11 = cos(qJ(3)) * pkin(3) + t13;
t9 = pkin(2) + t11;
t44 = t28 * t9 + pkin(1) + t33;
t20 = qJ(6) + t22;
t15 = sin(t20);
t16 = cos(t20);
t29 = cos(qJ(1));
t35 = t29 * t16;
t27 = sin(qJ(1));
t38 = t27 * t28;
t5 = t15 * t38 + t35;
t36 = t29 * t15;
t6 = -t16 * t38 + t36;
t43 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t27 * t16 - t28 * t36;
t8 = t27 * t15 + t28 * t35;
t42 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t18 = sin(t22);
t41 = pkin(5) * t18;
t12 = -pkin(4) * sin(t25) - t41;
t10 = sin(qJ(3)) * pkin(3) - t12;
t40 = pkin(7) + t10;
t37 = t28 * t29;
t32 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
t31 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t9;
t30 = -t31 * t26 + t39 * t28;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t44 * t27 + t40 * t29, t30 * t29, -t10 * t37 + t27 * t11 + t42, t12 * t37 + t27 * t13 + t42 (-t18 * t37 + t19 * t27) * pkin(5) + t42, t42; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t40 * t27 + t44 * t29, t30 * t27, -t10 * t38 - t29 * t11 + t43, t12 * t38 - t29 * t13 + t43 (-t18 * t38 - t19 * t29) * pkin(5) + t43, t43; 0, t31 * t28 + t33 (-t10 + t32) * t26 (t12 + t32) * t26 (t32 - t41) * t26, t32 * t26;];
Ja_transl  = t1;
