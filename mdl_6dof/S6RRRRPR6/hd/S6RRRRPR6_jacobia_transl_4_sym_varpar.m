% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:30
% EndTime: 2019-02-26 22:33:30
% DurationCPUTime: 0.10s
% Computational Cost: add. (93->28), mult. (124->43), div. (0->0), fcn. (136->8), ass. (0->27)
t17 = cos(qJ(2));
t14 = sin(qJ(2));
t28 = r_i_i_C(3) + pkin(9) + pkin(8);
t23 = t28 * t14;
t16 = cos(qJ(3));
t9 = pkin(3) * t16 + pkin(2);
t32 = t17 * t9 + pkin(1) + t23;
t12 = qJ(3) + qJ(4);
t10 = sin(t12);
t11 = cos(t12);
t18 = cos(qJ(1));
t15 = sin(qJ(1));
t27 = t15 * t17;
t5 = t10 * t27 + t11 * t18;
t6 = t10 * t18 - t11 * t27;
t31 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t26 = t17 * t18;
t7 = -t10 * t26 + t11 * t15;
t8 = t10 * t15 + t11 * t26;
t30 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t13 = sin(qJ(3));
t29 = pkin(3) * t13;
t24 = pkin(7) + t29;
t22 = -r_i_i_C(1) * t10 - r_i_i_C(2) * t11;
t21 = r_i_i_C(1) * t11 - r_i_i_C(2) * t10 + t9;
t20 = -t21 * t14 + t28 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t32 * t15 + t24 * t18, t20 * t18 (-t13 * t26 + t15 * t16) * pkin(3) + t30, t30, 0, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t24 * t15 + t32 * t18, t20 * t15 (-t13 * t27 - t16 * t18) * pkin(3) + t31, t31, 0, 0; 0, t21 * t17 + t23 (t22 - t29) * t14, t22 * t14, 0, 0;];
Ja_transl  = t1;
