% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:12
% EndTime: 2019-02-26 21:41:12
% DurationCPUTime: 0.16s
% Computational Cost: add. (222->44), mult. (306->63), div. (0->0), fcn. (365->10), ass. (0->35)
t18 = sin(qJ(2));
t22 = cos(qJ(2));
t17 = sin(qJ(4));
t35 = pkin(4) * t17 + qJ(3);
t21 = cos(qJ(4));
t38 = t21 * pkin(4) + pkin(2) + pkin(3);
t25 = t35 * t18 + t38 * t22;
t43 = pkin(1) + t25;
t36 = qJ(4) + pkin(10);
t33 = sin(t36);
t34 = cos(t36);
t8 = t18 * t34 - t22 * t33;
t16 = sin(qJ(6));
t20 = cos(qJ(6));
t27 = t20 * r_i_i_C(1) - t16 * r_i_i_C(2) + pkin(5);
t19 = sin(qJ(1));
t3 = t8 * t19;
t39 = -r_i_i_C(3) - pkin(9);
t7 = t18 * t33 + t22 * t34;
t4 = t7 * t19;
t42 = -t27 * t3 + t39 * t4;
t23 = cos(qJ(1));
t28 = t23 * t33;
t29 = t23 * t34;
t5 = -t18 * t28 - t22 * t29;
t6 = -t18 * t29 + t22 * t28;
t41 = t27 * t6 - t39 * t5;
t40 = -t27 * t7 - t39 * t8;
t37 = pkin(7) - qJ(5) - pkin(8);
t32 = -t16 * r_i_i_C(1) - t20 * r_i_i_C(2);
t26 = pkin(4) * (-t17 * t22 + t18 * t21);
t24 = -t38 * t18 + t35 * t22;
t2 = -t19 * t16 - t5 * t20;
t1 = t5 * t16 - t19 * t20;
t9 = [-t39 * t3 - t27 * t4 + (t32 + t37) * t23 - t43 * t19, t24 * t23 + t41, t23 * t18, t23 * t26 - t41, -t19, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); -t5 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t37 * t19 + t43 * t23 - t39 * t6, t24 * t19 + t42, t19 * t18, t19 * t26 - t42, t23 (-t4 * t16 + t23 * t20) * r_i_i_C(1) + (-t23 * t16 - t4 * t20) * r_i_i_C(2); 0, t25 - t40, -t22 (-t17 * t18 - t21 * t22) * pkin(4) + t40, 0, t32 * t8;];
Ja_transl  = t9;
