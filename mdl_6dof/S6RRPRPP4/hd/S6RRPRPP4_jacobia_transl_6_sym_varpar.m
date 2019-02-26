% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:36:35
% EndTime: 2019-02-26 21:36:36
% DurationCPUTime: 0.11s
% Computational Cost: add. (139->32), mult. (185->45), div. (0->0), fcn. (212->8), ass. (0->26)
t12 = sin(qJ(2));
t15 = cos(qJ(2));
t21 = pkin(2) + r_i_i_C(2) + qJ(5) + pkin(8);
t19 = t21 * t15;
t11 = sin(qJ(4));
t20 = pkin(4) * t11 + qJ(3);
t29 = -t20 * t12 - pkin(1) - t19;
t27 = pkin(5) + r_i_i_C(1);
t14 = cos(qJ(4));
t25 = pkin(4) * t14;
t26 = pkin(7) + pkin(3) + t25;
t13 = sin(qJ(1));
t24 = t13 * t12;
t16 = cos(qJ(1));
t23 = t16 * t12;
t22 = r_i_i_C(3) + qJ(6);
t9 = qJ(4) + pkin(9);
t7 = sin(t9);
t8 = cos(t9);
t18 = -t22 * t8 + t27 * t7 + t20;
t17 = -t21 * t12 + t18 * t15;
t4 = t16 * t8 - t7 * t24;
t3 = t16 * t7 + t8 * t24;
t2 = t13 * t8 + t7 * t23;
t1 = t13 * t7 - t8 * t23;
t5 = [t29 * t13 + t26 * t16 + t22 * t3 + t27 * t4, t17 * t16, t23, t22 * t2 - t27 * t1 + (-t11 * t13 + t14 * t23) * pkin(4), t16 * t15, t1; t22 * t1 + t26 * t13 - t29 * t16 + t27 * t2, t17 * t13, t24, -t22 * t4 + t27 * t3 + (t11 * t16 + t14 * t24) * pkin(4), t13 * t15, -t3; 0, t18 * t12 + t19, -t15 (-t22 * t7 - t27 * t8 - t25) * t15, t12, t15 * t8;];
Ja_transl  = t5;
