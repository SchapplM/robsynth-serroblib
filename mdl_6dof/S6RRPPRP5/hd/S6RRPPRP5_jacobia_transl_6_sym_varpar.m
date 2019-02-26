% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:27:26
% EndTime: 2019-02-26 21:27:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (134->28), mult. (173->39), div. (0->0), fcn. (200->8), ass. (0->23)
t12 = sin(qJ(2));
t14 = cos(qJ(2));
t20 = pkin(2) + r_i_i_C(2) + pkin(8) + qJ(4);
t18 = t20 * t14;
t19 = pkin(4) * sin(pkin(9)) + qJ(3);
t27 = -t19 * t12 - pkin(1) - t18;
t25 = pkin(5) + r_i_i_C(1);
t24 = pkin(7) + cos(pkin(9)) * pkin(4) + pkin(3);
t13 = sin(qJ(1));
t23 = t13 * t12;
t15 = cos(qJ(1));
t22 = t15 * t12;
t21 = r_i_i_C(3) + qJ(6);
t9 = pkin(9) + qJ(5);
t7 = sin(t9);
t8 = cos(t9);
t17 = -t21 * t8 + t25 * t7 + t19;
t16 = -t20 * t12 + t17 * t14;
t4 = t15 * t8 - t7 * t23;
t3 = t15 * t7 + t8 * t23;
t2 = t13 * t8 + t7 * t22;
t1 = t13 * t7 - t8 * t22;
t5 = [t27 * t13 + t24 * t15 + t21 * t3 + t25 * t4, t16 * t15, t22, t15 * t14, -t25 * t1 + t21 * t2, t1; t21 * t1 + t24 * t13 - t27 * t15 + t25 * t2, t16 * t13, t23, t13 * t14, -t21 * t4 + t25 * t3, -t3; 0, t17 * t12 + t18, -t14, t12 (-t21 * t7 - t25 * t8) * t14, t14 * t8;];
Ja_transl  = t5;
