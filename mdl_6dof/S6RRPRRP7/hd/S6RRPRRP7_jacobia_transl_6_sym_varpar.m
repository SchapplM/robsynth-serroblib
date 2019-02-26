% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP7
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
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:25
% EndTime: 2019-02-26 21:49:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (168->36), mult. (414->52), div. (0->0), fcn. (509->8), ass. (0->31)
t19 = sin(qJ(5));
t21 = cos(qJ(5));
t34 = r_i_i_C(3) + qJ(6);
t40 = pkin(5) + r_i_i_C(1);
t25 = t34 * t19 + t40 * t21 + pkin(4);
t38 = pkin(9) + r_i_i_C(2);
t22 = cos(qJ(2));
t35 = sin(qJ(4));
t36 = sin(qJ(2));
t37 = cos(qJ(4));
t12 = -t22 * t35 + t36 * t37;
t20 = sin(qJ(1));
t7 = t12 * t20;
t11 = t22 * t37 + t36 * t35;
t8 = t11 * t20;
t48 = t25 * t7 + t38 * t8;
t23 = cos(qJ(1));
t10 = t12 * t23;
t9 = t11 * t23;
t47 = t25 * t10 + t38 * t9;
t46 = -t25 * t11 + t38 * t12;
t41 = pkin(2) + pkin(3);
t27 = t36 * qJ(3) + t41 * t22;
t43 = pkin(1) + t27;
t39 = pkin(7) - pkin(8);
t1 = t8 * t19 - t23 * t21;
t28 = t23 * t19 + t8 * t21;
t24 = qJ(3) * t22 - t41 * t36;
t6 = -t20 * t19 + t9 * t21;
t5 = t9 * t19 + t20 * t21;
t2 = [-t8 * pkin(4) - t34 * t1 - t43 * t20 + t39 * t23 - t40 * t28 + t38 * t7, t24 * t23 - t47, t23 * t36, t47, t34 * t6 - t40 * t5, t5; t9 * pkin(4) - t38 * t10 + t39 * t20 + t43 * t23 + t34 * t5 + t40 * t6, t24 * t20 - t48, t20 * t36, t48, -t40 * t1 + t34 * t28, t1; 0, t27 - t46, -t22, t46 (-t40 * t19 + t34 * t21) * t12, t12 * t19;];
Ja_transl  = t2;
