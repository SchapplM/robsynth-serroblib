% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR12_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR12_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:00:24
% EndTime: 2019-02-26 22:00:24
% DurationCPUTime: 0.14s
% Computational Cost: add. (150->40), mult. (261->61), div. (0->0), fcn. (320->10), ass. (0->33)
t21 = cos(pkin(6));
t23 = sin(qJ(2));
t27 = cos(qJ(1));
t33 = t27 * t23;
t24 = sin(qJ(1));
t26 = cos(qJ(2));
t34 = t24 * t26;
t11 = t21 * t34 + t33;
t19 = qJ(4) + qJ(5);
t17 = sin(t19);
t18 = cos(t19);
t20 = sin(pkin(6));
t38 = t20 * t24;
t5 = t11 * t18 - t17 * t38;
t6 = t11 * t17 + t18 * t38;
t42 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t36 = t20 * t27;
t13 = t18 * t36;
t32 = t27 * t26;
t35 = t24 * t23;
t9 = -t21 * t32 + t35;
t41 = (t17 * t36 + t9 * t18) * r_i_i_C(1) + (-t9 * t17 + t13) * r_i_i_C(2);
t37 = t20 * t26;
t40 = (-t21 * t17 - t18 * t37) * r_i_i_C(1) + (t17 * t37 - t21 * t18) * r_i_i_C(2);
t25 = cos(qJ(4));
t39 = t25 * pkin(4) + pkin(3) + pkin(8);
t31 = -r_i_i_C(3) - pkin(10) - pkin(9) - pkin(2);
t22 = sin(qJ(4));
t30 = t22 * pkin(4) + qJ(3);
t29 = t17 * r_i_i_C(1) + t18 * r_i_i_C(2) + t30;
t12 = -t21 * t35 + t32;
t10 = t21 * t33 + t34;
t1 = [-t24 * pkin(1) + t13 * r_i_i_C(1) - t29 * t9 + (-r_i_i_C(2) * t17 + t39) * t36 + t31 * t10, t11 * t31 + t12 * t29, t11 (t11 * t25 - t22 * t38) * pkin(4) + t42, t42, 0; t27 * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t11 * t30 - t12 * t31 + t38 * t39, t10 * t29 + t31 * t9, t9 (t22 * t36 + t25 * t9) * pkin(4) + t41, t41, 0; 0 (t23 * t29 - t26 * t31) * t20, -t37 (-t21 * t22 - t25 * t37) * pkin(4) + t40, t40, 0;];
Ja_transl  = t1;
