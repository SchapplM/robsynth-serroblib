% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:47
% EndTime: 2019-02-26 21:03:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (205->35), mult. (140->46), div. (0->0), fcn. (157->11), ass. (0->29)
t21 = pkin(10) + qJ(3);
t17 = cos(t21);
t16 = sin(t21);
t35 = r_i_i_C(3) + pkin(9) + qJ(5) + pkin(8);
t30 = t35 * t16;
t22 = qJ(4) + pkin(11);
t11 = pkin(5) * cos(t22) + cos(qJ(4)) * pkin(4);
t9 = pkin(3) + t11;
t38 = t30 + t17 * t9 + cos(pkin(10)) * pkin(2) + pkin(1);
t18 = qJ(6) + t22;
t14 = cos(t18);
t25 = cos(qJ(1));
t13 = sin(t18);
t24 = sin(qJ(1));
t33 = t24 * t13;
t5 = t14 * t25 + t17 * t33;
t32 = t24 * t14;
t6 = t13 * t25 - t17 * t32;
t37 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t34 = t17 * t25;
t7 = -t13 * t34 + t32;
t8 = t14 * t34 + t33;
t36 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t10 = pkin(5) * sin(t22) + sin(qJ(4)) * pkin(4);
t31 = t10 + pkin(7) + qJ(2);
t28 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t27 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t26 = -t27 * t16 + t35 * t17;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t38 * t24 + t31 * t25, t24, t26 * t25, -t10 * t34 + t24 * t11 + t36, t25 * t16, t36; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t31 * t24 + t38 * t25, -t25, t26 * t24, -t24 * t17 * t10 - t11 * t25 + t37, t24 * t16, t37; 0, 0, t27 * t17 + t30 (-t10 + t28) * t16, -t17, t28 * t16;];
Ja_transl  = t1;
