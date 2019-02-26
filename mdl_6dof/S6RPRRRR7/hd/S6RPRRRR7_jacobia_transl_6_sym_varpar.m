% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:18:10
% EndTime: 2019-02-26 21:18:10
% DurationCPUTime: 0.14s
% Computational Cost: add. (209->41), mult. (166->46), div. (0->0), fcn. (176->10), ass. (0->34)
t21 = qJ(3) + qJ(4);
t19 = qJ(5) + t21;
t15 = sin(t19);
t45 = pkin(10) + r_i_i_C(3);
t47 = t45 * t15;
t16 = cos(t19);
t46 = t45 * t16;
t43 = sin(qJ(3)) * pkin(3);
t42 = pkin(4) * sin(t21);
t41 = pkin(4) * cos(t21);
t22 = sin(qJ(6));
t40 = r_i_i_C(2) * t22;
t39 = pkin(1) + pkin(9) + pkin(8) + pkin(7);
t26 = cos(qJ(1));
t37 = t22 * t26;
t24 = sin(qJ(1));
t36 = t24 * t22;
t25 = cos(qJ(6));
t35 = t24 * t25;
t34 = t25 * t26;
t33 = t16 * t40;
t32 = -r_i_i_C(1) * t25 - pkin(5);
t31 = t24 * t47 + (pkin(5) * t24 + r_i_i_C(1) * t35) * t16;
t30 = t46 + (t32 + t40) * t15;
t29 = pkin(5) * t15 + qJ(2) + t42 + t43 - t46;
t28 = t32 * t16 - t47;
t27 = t30 - t42;
t8 = t41 + cos(qJ(3)) * pkin(3);
t6 = t26 * t33;
t4 = t15 * t34 - t36;
t3 = t15 * t37 + t35;
t2 = t15 * t35 + t37;
t1 = -t15 * t36 + t34;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t39 * t24 + t29 * t26, t24 (t8 - t33) * t24 + t31 (-t33 + t41) * t24 + t31, -t24 * t33 + t31, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t29 * t24 + t39 * t26, -t26, t6 + (-t8 + t28) * t26, t6 + (t28 - t41) * t26, t28 * t26 + t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t27 - t43, t27, t30 (-r_i_i_C(1) * t22 - r_i_i_C(2) * t25) * t16;];
Ja_transl  = t5;
