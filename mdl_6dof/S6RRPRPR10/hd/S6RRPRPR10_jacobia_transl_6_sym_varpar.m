% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR10_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR10_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:00
% EndTime: 2019-02-26 21:43:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (294->49), mult. (477->80), div. (0->0), fcn. (606->12), ass. (0->36)
t18 = cos(pkin(11)) * pkin(3) + pkin(2);
t21 = pkin(11) + qJ(4);
t19 = sin(t21);
t20 = cos(t21);
t25 = sin(qJ(6));
t28 = cos(qJ(6));
t34 = t25 * r_i_i_C(1) + t28 * r_i_i_C(2) + qJ(5);
t38 = pkin(4) + pkin(10) + r_i_i_C(3);
t45 = t34 * t19 + t38 * t20 + t18;
t44 = pkin(5) + pkin(9) + qJ(3);
t23 = sin(pkin(6));
t26 = sin(qJ(2));
t43 = t23 * t26;
t27 = sin(qJ(1));
t42 = t23 * t27;
t29 = cos(qJ(2));
t41 = t23 * t29;
t30 = cos(qJ(1));
t40 = t23 * t30;
t39 = cos(pkin(6));
t37 = t27 * t39;
t36 = t30 * t39;
t35 = t23 * (pkin(3) * sin(pkin(11)) + pkin(8));
t12 = t26 * t36 + t27 * t29;
t3 = t12 * t19 + t20 * t40;
t33 = -t12 * t20 + t19 * t40;
t32 = t28 * r_i_i_C(1) - t25 * r_i_i_C(2) + t44;
t14 = -t26 * t37 + t30 * t29;
t13 = t30 * t26 + t29 * t37;
t11 = t27 * t26 - t29 * t36;
t9 = t19 * t43 - t39 * t20;
t8 = t14 * t20 + t19 * t42;
t7 = t14 * t19 - t20 * t42;
t2 = t13 * t28 + t7 * t25;
t1 = -t13 * t25 + t7 * t28;
t4 = [-t27 * pkin(1) - t32 * t11 - t12 * t18 - t34 * t3 + t30 * t35 + t38 * t33, -t13 * t45 + t32 * t14, t13, t34 * t8 - t38 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t30 * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t7 * qJ(5) + t44 * t13 + t14 * t18 + t27 * t35 + t38 * t8, -t11 * t45 + t32 * t12, t11, -t38 * t3 - t34 * t33, t3 (-t11 * t25 + t3 * t28) * r_i_i_C(1) + (-t11 * t28 - t3 * t25) * r_i_i_C(2); 0 (t32 * t26 + t45 * t29) * t23, -t41, -t38 * t9 + t34 * (t39 * t19 + t20 * t43) t9 (t25 * t41 + t9 * t28) * r_i_i_C(1) + (-t9 * t25 + t28 * t41) * r_i_i_C(2);];
Ja_transl  = t4;
