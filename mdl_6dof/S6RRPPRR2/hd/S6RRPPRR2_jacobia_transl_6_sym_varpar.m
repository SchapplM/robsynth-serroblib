% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:28:44
% EndTime: 2019-02-26 21:28:44
% DurationCPUTime: 0.15s
% Computational Cost: add. (203->37), mult. (140->49), div. (0->0), fcn. (157->12), ass. (0->30)
t22 = qJ(2) + pkin(10);
t15 = sin(t22);
t36 = r_i_i_C(3) + pkin(9) + pkin(8) + qJ(4);
t41 = cos(qJ(2)) * pkin(2) + t36 * t15;
t17 = cos(t22);
t21 = pkin(11) + qJ(5);
t16 = cos(t21);
t9 = pkin(5) * t16 + cos(pkin(11)) * pkin(4) + pkin(3);
t40 = t17 * t9 + pkin(1) + t41;
t18 = qJ(6) + t21;
t12 = cos(t18);
t26 = cos(qJ(1));
t11 = sin(t18);
t25 = sin(qJ(1));
t34 = t25 * t11;
t5 = t12 * t26 + t17 * t34;
t33 = t25 * t12;
t6 = t11 * t26 - t17 * t33;
t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t35 = t17 * t26;
t7 = -t11 * t35 + t33;
t8 = t12 * t35 + t34;
t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t14 = sin(t21);
t37 = pkin(5) * t14;
t32 = t37 + sin(pkin(11)) * pkin(4) + qJ(3) + pkin(7);
t29 = -r_i_i_C(1) * t11 - r_i_i_C(2) * t12;
t28 = r_i_i_C(1) * t12 - r_i_i_C(2) * t11 + t9;
t27 = -sin(qJ(2)) * pkin(2) - t15 * t28 + t17 * t36;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t25 + t32 * t26, t27 * t26, t25, t26 * t15 (-t14 * t35 + t16 * t25) * pkin(5) + t38, t38; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t32 * t25 + t40 * t26, t27 * t25, -t26, t25 * t15 (-t14 * t17 * t25 - t16 * t26) * pkin(5) + t39, t39; 0, t17 * t28 + t41, 0, -t17 (t29 - t37) * t15, t29 * t15;];
Ja_transl  = t1;
