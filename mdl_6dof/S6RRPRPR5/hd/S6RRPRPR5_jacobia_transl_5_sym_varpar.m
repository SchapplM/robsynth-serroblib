% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_jacobia_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:07
% EndTime: 2019-02-26 21:40:08
% DurationCPUTime: 0.19s
% Computational Cost: add. (255->52), mult. (641->83), div. (0->0), fcn. (844->12), ass. (0->37)
t38 = cos(qJ(2));
t51 = t38 * pkin(2);
t33 = cos(pkin(6));
t50 = t33 * t38;
t30 = sin(pkin(6));
t36 = sin(qJ(1));
t49 = t36 * t30;
t39 = cos(qJ(1));
t48 = t39 * t30;
t47 = r_i_i_C(3) + qJ(5);
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t29 = sin(pkin(11));
t32 = cos(pkin(11));
t35 = sin(qJ(2));
t44 = t38 * t29 + t35 * t32;
t19 = t44 * t33;
t21 = t35 * t29 - t38 * t32;
t9 = t39 * t19 - t36 * t21;
t46 = -t34 * t48 + t9 * t37;
t45 = t36 * t19 + t39 * t21;
t28 = sin(pkin(12));
t31 = cos(pkin(12));
t43 = r_i_i_C(1) * t31 - r_i_i_C(2) * t28 + pkin(4);
t42 = t28 * r_i_i_C(1) + t31 * r_i_i_C(2) + pkin(9);
t1 = t9 * t34 + t37 * t48;
t41 = t21 * t33;
t40 = t47 * t34 + t43 * t37 + pkin(3);
t27 = pkin(1) + t51;
t20 = t33 * t35 * pkin(2) + (-pkin(8) - qJ(3)) * t30;
t18 = t44 * t30;
t13 = t18 * t34 - t33 * t37;
t11 = t36 * t41 - t39 * t44;
t8 = -t36 * t44 - t39 * t41;
t6 = t34 * t49 - t37 * t45;
t5 = -t34 * t45 - t37 * t49;
t2 = [(t8 * t28 - t31 * t46) * r_i_i_C(1) + (t28 * t46 + t8 * t31) * r_i_i_C(2) - t46 * pkin(4) - t9 * pkin(3) + t8 * pkin(9) - t36 * t27 - t39 * t20 - t47 * t1 (-t39 * t35 - t36 * t50) * pkin(2) - t42 * t45 + t40 * t11, t49, -t43 * t5 + t47 * t6, t5, 0; (-t11 * t28 + t6 * t31) * r_i_i_C(1) + (-t11 * t31 - t6 * t28) * r_i_i_C(2) + t6 * pkin(4) - t45 * pkin(3) - t11 * pkin(9) + t39 * t27 - t36 * t20 + t47 * t5 (-t36 * t35 + t39 * t50) * pkin(2) + t42 * t9 + t40 * t8, -t48, -t43 * t1 + t47 * t46, t1, 0; 0, t42 * t18 + (-t21 * t40 + t51) * t30, t33, t47 * (t18 * t37 + t33 * t34) - t43 * t13, t13, 0;];
Ja_transl  = t2;
