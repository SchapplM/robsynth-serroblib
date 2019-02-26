% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:47:14
% EndTime: 2019-02-26 21:47:14
% DurationCPUTime: 0.13s
% Computational Cost: add. (141->32), mult. (131->44), div. (0->0), fcn. (145->10), ass. (0->31)
t16 = qJ(2) + pkin(10);
t11 = sin(t16);
t36 = r_i_i_C(3) + pkin(9) + pkin(8);
t41 = cos(qJ(2)) * pkin(2) + t36 * t11;
t12 = cos(t16);
t22 = cos(qJ(4));
t9 = pkin(4) * t22 + pkin(3);
t40 = t12 * t9 + pkin(1) + t41;
t17 = qJ(4) + qJ(5);
t14 = cos(t17);
t23 = cos(qJ(1));
t31 = t14 * t23;
t13 = sin(t17);
t21 = sin(qJ(1));
t34 = t13 * t21;
t5 = t12 * t34 + t31;
t32 = t14 * t21;
t33 = t13 * t23;
t6 = -t12 * t32 + t33;
t39 = -t5 * r_i_i_C(1) + t6 * r_i_i_C(2);
t7 = -t12 * t33 + t32;
t8 = t12 * t31 + t34;
t38 = t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t19 = sin(qJ(4));
t37 = pkin(4) * t19;
t35 = t12 * t19;
t29 = qJ(3) + pkin(7) + t37;
t27 = -r_i_i_C(1) * t13 - r_i_i_C(2) * t14;
t26 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + t9;
t25 = -sin(qJ(2)) * pkin(2) - t11 * t26 + t12 * t36;
t1 = [t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t40 * t21 + t29 * t23, t25 * t23, t21 (t21 * t22 - t23 * t35) * pkin(4) + t38, t38, 0; t8 * r_i_i_C(1) + t7 * r_i_i_C(2) + t29 * t21 + t40 * t23, t25 * t21, -t23 (-t21 * t35 - t22 * t23) * pkin(4) + t39, t39, 0; 0, t12 * t26 + t41, 0 (t27 - t37) * t11, t27 * t11, 0;];
Ja_transl  = t1;
