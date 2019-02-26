% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% Ja_transl [3x7]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S7RRRRRRR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobia_transl_4_sym_varpar: qJ has to be [7x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S7RRRRRRR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobia_transl_4_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:20
% EndTime: 2019-02-26 22:54:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (56->31), mult. (157->63), div. (0->0), fcn. (189->8), ass. (0->29)
t14 = cos(qJ(2));
t9 = sin(qJ(3));
t29 = t14 * t9;
t15 = cos(qJ(1));
t28 = t15 * t9;
t10 = sin(qJ(2));
t11 = sin(qJ(1));
t27 = t10 * t11;
t12 = cos(qJ(4));
t26 = t10 * t12;
t13 = cos(qJ(3));
t25 = t10 * t13;
t24 = t10 * t15;
t23 = t13 * t14;
t22 = t15 * t13;
t8 = sin(qJ(4));
t21 = r_i_i_C(1) * t12 - r_i_i_C(2) * t8;
t4 = t11 * t23 + t28;
t20 = -t12 * t4 - t8 * t27;
t19 = t11 * t26 - t4 * t8;
t18 = t12 * t14 + t8 * t25;
t17 = -t12 * t25 + t14 * t8;
t16 = t10 * t9 * r_i_i_C(3) - t14 * pkin(2) + t17 * r_i_i_C(1) + t18 * r_i_i_C(2);
t6 = -t11 * t9 + t14 * t22;
t5 = -t11 * t13 - t14 * t28;
t3 = t11 * t29 - t22;
t2 = t6 * t12 + t8 * t24;
t1 = t12 * t24 - t6 * t8;
t7 = [pkin(2) * t27 + t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t3 * r_i_i_C(3), t16 * t15, -r_i_i_C(3) * t6 + t21 * t5, r_i_i_C(1) * t1 - t2 * r_i_i_C(2), 0, 0, 0; -pkin(2) * t24 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t5 * r_i_i_C(3), t16 * t11, -r_i_i_C(3) * t4 - t21 * t3, t19 * r_i_i_C(1) + t20 * r_i_i_C(2), 0, 0, 0; 0 (t10 * t8 + t12 * t23) * r_i_i_C(1) + (-t8 * t23 + t26) * r_i_i_C(2) - r_i_i_C(3) * t29 - t10 * pkin(2) (-r_i_i_C(3) * t13 - t21 * t9) * t10, -t18 * r_i_i_C(1) + t17 * r_i_i_C(2), 0, 0, 0;];
Ja_transl  = t7;
