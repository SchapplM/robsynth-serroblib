% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function Ja_transl = S6RRRRPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:33:30
% EndTime: 2019-02-26 22:33:30
% DurationCPUTime: 0.13s
% Computational Cost: add. (294->40), mult. (187->52), div. (0->0), fcn. (206->12), ass. (0->32)
t27 = cos(qJ(2));
t25 = sin(qJ(2));
t38 = r_i_i_C(3) + pkin(10) + qJ(5) + pkin(9) + pkin(8);
t32 = t38 * t25;
t24 = qJ(3) + qJ(4);
t20 = pkin(11) + t24;
t13 = pkin(5) * cos(t20) + pkin(4) * cos(t24);
t11 = cos(qJ(3)) * pkin(3) + t13;
t9 = pkin(2) + t11;
t42 = t27 * t9 + pkin(1) + t32;
t19 = qJ(6) + t20;
t15 = sin(t19);
t16 = cos(t19);
t28 = cos(qJ(1));
t34 = t28 * t16;
t26 = sin(qJ(1));
t37 = t26 * t27;
t5 = t15 * t37 + t34;
t35 = t28 * t15;
t6 = -t16 * t37 + t35;
t41 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t26 * t16 - t27 * t35;
t8 = t26 * t15 + t27 * t34;
t40 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t12 = -pkin(4) * sin(t24) - pkin(5) * sin(t20);
t10 = sin(qJ(3)) * pkin(3) - t12;
t39 = pkin(7) + t10;
t36 = t27 * t28;
t31 = -r_i_i_C(1) * t15 - r_i_i_C(2) * t16;
t30 = r_i_i_C(1) * t16 - r_i_i_C(2) * t15 + t9;
t29 = -t30 * t25 + t38 * t27;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t42 * t26 + t39 * t28, t29 * t28, -t10 * t36 + t26 * t11 + t40, t12 * t36 + t26 * t13 + t40, t28 * t25, t40; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t39 * t26 + t42 * t28, t29 * t26, -t10 * t37 - t28 * t11 + t41, t12 * t37 - t28 * t13 + t41, t26 * t25, t41; 0, t30 * t27 + t32 (-t10 + t31) * t25 (t12 + t31) * t25, -t27, t31 * t25;];
Ja_transl  = t1;
