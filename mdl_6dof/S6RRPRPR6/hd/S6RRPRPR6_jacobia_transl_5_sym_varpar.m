% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (204->42), mult. (503->67), div. (0->0), fcn. (658->10), ass. (0->35)
t45 = pkin(4) - r_i_i_C(2);
t44 = pkin(9) + r_i_i_C(1);
t33 = cos(qJ(2));
t43 = t33 * pkin(2);
t28 = cos(pkin(6));
t42 = t28 * t33;
t26 = sin(pkin(6));
t31 = sin(qJ(1));
t41 = t31 * t26;
t34 = cos(qJ(1));
t40 = t34 * t26;
t39 = r_i_i_C(3) + qJ(5);
t25 = sin(pkin(11));
t27 = cos(pkin(11));
t30 = sin(qJ(2));
t37 = t33 * t25 + t30 * t27;
t18 = t37 * t28;
t20 = t30 * t25 - t33 * t27;
t9 = t34 * t18 - t31 * t20;
t38 = t31 * t18 + t34 * t20;
t29 = sin(qJ(4));
t32 = cos(qJ(4));
t1 = t9 * t29 + t32 * t40;
t36 = t29 * t40 - t9 * t32;
t35 = t39 * t29 + t45 * t32 + pkin(3);
t24 = pkin(1) + t43;
t19 = t28 * t30 * pkin(2) + (-pkin(8) - qJ(3)) * t26;
t17 = t20 * t28;
t16 = t37 * t26;
t13 = t16 * t29 - t28 * t32;
t11 = t31 * t17 - t34 * t37;
t8 = -t34 * t17 - t31 * t37;
t6 = t29 * t41 - t32 * t38;
t5 = -t29 * t38 - t32 * t41;
t2 = [-t9 * pkin(3) - t39 * t1 - t34 * t19 - t31 * t24 + t45 * t36 + t44 * t8, -t44 * t38 + (-t30 * t34 - t31 * t42) * pkin(2) + t35 * t11, t41, t39 * t6 - t45 * t5, t5, 0; -pkin(3) * t38 - t44 * t11 - t31 * t19 + t34 * t24 + t39 * t5 + t45 * t6, t44 * t9 + (-t30 * t31 + t34 * t42) * pkin(2) + t35 * t8, -t40, -t45 * t1 - t39 * t36, t1, 0; 0, t44 * t16 + (-t20 * t35 + t43) * t26, t28, t39 * (t16 * t32 + t28 * t29) - t45 * t13, t13, 0;];
Ja_transl  = t2;
