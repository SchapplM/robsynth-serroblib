% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:57:34
% EndTime: 2019-02-26 21:57:35
% DurationCPUTime: 0.17s
% Computational Cost: add. (190->43), mult. (339->62), div. (0->0), fcn. (409->10), ass. (0->38)
t50 = pkin(7) - pkin(8);
t23 = sin(qJ(2));
t26 = cos(qJ(2));
t45 = pkin(2) + pkin(3);
t30 = t23 * qJ(3) + t45 * t26;
t49 = pkin(1) + t30;
t39 = sin(qJ(4));
t40 = cos(qJ(4));
t12 = t23 * t40 - t26 * t39;
t25 = cos(qJ(5));
t18 = t25 * pkin(5) + pkin(4);
t21 = qJ(5) + qJ(6);
t19 = sin(t21);
t20 = cos(t21);
t32 = t20 * r_i_i_C(1) - t19 * r_i_i_C(2) + t18;
t41 = r_i_i_C(3) + pkin(10) + pkin(9);
t24 = sin(qJ(1));
t7 = t12 * t24;
t11 = t23 * t39 + t26 * t40;
t8 = t11 * t24;
t48 = t32 * t7 + t41 * t8;
t27 = cos(qJ(1));
t34 = t27 * t39;
t35 = t27 * t40;
t10 = -t23 * t35 + t26 * t34;
t9 = -t23 * t34 - t26 * t35;
t47 = -t32 * t10 - t41 * t9;
t46 = -t32 * t11 + t41 * t12;
t44 = (-t8 * t19 + t27 * t20) * r_i_i_C(1) + (-t27 * t19 - t8 * t20) * r_i_i_C(2);
t5 = t9 * t19 - t24 * t20;
t6 = -t24 * t19 - t9 * t20;
t43 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t22 = sin(qJ(5));
t42 = t22 * pkin(5);
t33 = -t19 * r_i_i_C(1) - t20 * r_i_i_C(2);
t31 = qJ(3) * t26 - t45 * t23;
t29 = t33 - t42;
t1 = [t41 * t7 - t32 * t8 + (t29 + t50) * t27 - t49 * t24, t31 * t27 - t47, t27 * t23, t47 (t22 * t9 - t24 * t25) * pkin(5) + t43, t43; t6 * r_i_i_C(1) + t5 * r_i_i_C(2) - t9 * t18 + t41 * t10 + (-t42 + t50) * t24 + t49 * t27, t31 * t24 - t48, t24 * t23, t48 (-t22 * t8 + t25 * t27) * pkin(5) + t44, t44; 0, t30 - t46, -t26, t46, t29 * t12, t33 * t12;];
Ja_transl  = t1;
