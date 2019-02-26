% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.14s
% Computational Cost: add. (216->39), mult. (290->55), div. (0->0), fcn. (349->10), ass. (0->32)
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t33 = pkin(4) * sin(pkin(10)) + qJ(3);
t36 = pkin(2) + cos(pkin(10)) * pkin(4) + pkin(3);
t24 = t33 * t18 + t36 * t21;
t41 = pkin(1) + t24;
t34 = pkin(10) + qJ(5);
t31 = sin(t34);
t32 = cos(t34);
t8 = t18 * t32 - t21 * t31;
t17 = sin(qJ(6));
t20 = cos(qJ(6));
t25 = t20 * r_i_i_C(1) - t17 * r_i_i_C(2) + pkin(5);
t19 = sin(qJ(1));
t3 = t8 * t19;
t37 = -r_i_i_C(3) - pkin(9);
t7 = t18 * t31 + t21 * t32;
t4 = t7 * t19;
t40 = t25 * t3 - t37 * t4;
t22 = cos(qJ(1));
t26 = t22 * t31;
t27 = t22 * t32;
t5 = -t18 * t26 - t21 * t27;
t6 = -t18 * t27 + t21 * t26;
t39 = -t25 * t6 + t37 * t5;
t38 = -t25 * t7 - t37 * t8;
t35 = pkin(7) - pkin(8) - qJ(4);
t30 = -t17 * r_i_i_C(1) - t20 * r_i_i_C(2);
t23 = -t36 * t18 + t33 * t21;
t2 = -t19 * t17 - t5 * t20;
t1 = t5 * t17 - t19 * t20;
t9 = [-t37 * t3 - t25 * t4 + (t30 + t35) * t22 - t41 * t19, t23 * t22 - t39, t22 * t18, -t19, t39, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t35 * t19 + t41 * t22 - t37 * t6, t23 * t19 - t40, t19 * t18, t22, t40 (-t4 * t17 + t22 * t20) * r_i_i_C(1) + (-t22 * t17 - t4 * t20) * r_i_i_C(2); 0, t24 - t38, -t21, 0, t38, t30 * t8;];
Ja_transl  = t9;
