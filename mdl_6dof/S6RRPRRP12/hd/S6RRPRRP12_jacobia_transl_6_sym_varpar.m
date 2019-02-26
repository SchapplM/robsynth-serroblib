% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP12_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP12_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:13
% EndTime: 2019-02-26 21:52:13
% DurationCPUTime: 0.13s
% Computational Cost: add. (176->32), mult. (223->44), div. (0->0), fcn. (255->8), ass. (0->29)
t36 = pkin(5) + r_i_i_C(1);
t31 = r_i_i_C(3) + qJ(6);
t17 = sin(qJ(2));
t20 = cos(qJ(2));
t30 = pkin(2) + r_i_i_C(2) + pkin(9) + pkin(8);
t28 = t30 * t20;
t16 = sin(qJ(4));
t29 = pkin(4) * t16 + qJ(3);
t38 = -t29 * t17 - pkin(1) - t28;
t19 = cos(qJ(4));
t35 = t19 * pkin(4);
t34 = pkin(7) + pkin(3) + t35;
t18 = sin(qJ(1));
t33 = t18 * t17;
t21 = cos(qJ(1));
t32 = t21 * t17;
t15 = qJ(4) + qJ(5);
t13 = sin(t15);
t14 = cos(t15);
t7 = t18 * t13 - t14 * t32;
t8 = t13 * t32 + t18 * t14;
t27 = t31 * t8 - t36 * t7;
t10 = -t13 * t33 + t21 * t14;
t9 = t21 * t13 + t14 * t33;
t26 = -t31 * t10 + t36 * t9;
t25 = -t31 * t13 - t36 * t14;
t24 = t36 * t13 - t31 * t14 + t29;
t23 = -t30 * t17 + t24 * t20;
t1 = [t36 * t10 + t38 * t18 + t34 * t21 + t31 * t9, t23 * t21, t32 (-t16 * t18 + t19 * t32) * pkin(4) + t27, t27, t7; t34 * t18 - t38 * t21 + t31 * t7 + t36 * t8, t23 * t18, t33 (t16 * t21 + t19 * t33) * pkin(4) + t26, t26, -t9; 0, t24 * t17 + t28, -t20 (t25 - t35) * t20, t25 * t20, t20 * t14;];
Ja_transl  = t1;
