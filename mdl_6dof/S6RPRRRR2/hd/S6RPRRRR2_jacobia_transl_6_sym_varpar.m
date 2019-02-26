% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR2_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR2_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:15:25
% EndTime: 2019-02-26 21:15:25
% DurationCPUTime: 0.16s
% Computational Cost: add. (243->40), mult. (167->55), div. (0->0), fcn. (179->12), ass. (0->35)
t26 = qJ(5) + qJ(6);
t20 = sin(t26);
t27 = qJ(3) + qJ(4);
t21 = sin(t27);
t23 = cos(t27);
t52 = r_i_i_C(2) * t20 * t21 + r_i_i_C(3) * t23;
t30 = cos(qJ(5));
t16 = t30 * pkin(5) + pkin(4);
t31 = -pkin(10) - pkin(9);
t51 = t23 * t16 + (r_i_i_C(3) - t31) * t21;
t24 = cos(qJ(3)) * pkin(3);
t50 = pkin(2) + t24 + t51;
t25 = qJ(1) + pkin(11);
t18 = sin(t25);
t19 = cos(t25);
t22 = cos(t26);
t43 = t20 * t23;
t5 = t18 * t43 + t19 * t22;
t42 = t22 * t23;
t6 = -t18 * t42 + t19 * t20;
t49 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = t18 * t22 - t19 * t43;
t8 = t18 * t20 + t19 * t42;
t48 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t28 = sin(qJ(5));
t47 = pkin(5) * t28;
t44 = t52 * t18;
t41 = t23 * t28;
t40 = t52 * t19;
t38 = pkin(8) + pkin(7) + t47;
t36 = -r_i_i_C(1) * t20 - r_i_i_C(2) * t22;
t35 = -t23 * t31 + (-r_i_i_C(1) * t22 - t16) * t21;
t34 = r_i_i_C(1) * t42 - r_i_i_C(2) * t43 + t51;
t33 = -sin(qJ(3)) * pkin(3) + t35;
t1 = [-sin(qJ(1)) * pkin(1) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t38 * t19 - t50 * t18, 0, t19 * t33 + t40, t19 * t35 + t40 (t18 * t30 - t19 * t41) * pkin(5) + t48, t48; cos(qJ(1)) * pkin(1) + t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t38 * t18 + t50 * t19, 0, t18 * t33 + t44, t18 * t35 + t44 (-t18 * t41 - t19 * t30) * pkin(5) + t49, t49; 0, 1, t24 + t34, t34 (t36 - t47) * t21, t36 * t21;];
Ja_transl  = t1;
