% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPP1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:51
% EndTime: 2019-02-26 21:34:51
% DurationCPUTime: 0.12s
% Computational Cost: add. (122->33), mult. (113->46), div. (0->0), fcn. (126->10), ass. (0->25)
t28 = r_i_i_C(3) + qJ(5) + pkin(8);
t13 = qJ(2) + pkin(9);
t8 = sin(t13);
t31 = cos(qJ(2)) * pkin(2) + t28 * t8;
t10 = cos(t13);
t19 = cos(qJ(4));
t5 = pkin(4) * t19 + pkin(3);
t30 = t10 * t5 + pkin(1) + t31;
t16 = sin(qJ(4));
t29 = pkin(4) * t16;
t18 = sin(qJ(1));
t27 = t10 * t18;
t20 = cos(qJ(1));
t26 = t10 * t20;
t23 = qJ(3) + pkin(7) + t29;
t12 = qJ(4) + pkin(10);
t7 = sin(t12);
t9 = cos(t12);
t22 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + t5;
t21 = -sin(qJ(2)) * pkin(2) + t28 * t10 - t22 * t8;
t4 = t18 * t7 + t9 * t26;
t3 = t18 * t9 - t7 * t26;
t2 = t20 * t7 - t9 * t27;
t1 = t20 * t9 + t7 * t27;
t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t30 * t18 + t23 * t20, t21 * t20, t18, t3 * r_i_i_C(1) - t4 * r_i_i_C(2) + (-t16 * t26 + t18 * t19) * pkin(4), t20 * t8, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t23 * t18 + t30 * t20, t21 * t18, -t20, -t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + (-t16 * t27 - t19 * t20) * pkin(4), t18 * t8, 0; 0, t22 * t10 + t31, 0 (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9 - t29) * t8, -t10, 0;];
Ja_transl  = t6;
