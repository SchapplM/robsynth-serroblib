% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:17:17
% EndTime: 2019-02-26 22:17:17
% DurationCPUTime: 0.13s
% Computational Cost: add. (144->27), mult. (157->37), div. (0->0), fcn. (182->8), ass. (0->30)
t23 = qJ(2) + qJ(3);
t20 = sin(t23);
t21 = cos(t23);
t41 = pkin(3) + pkin(4);
t46 = t20 * qJ(4) + t41 * t21;
t27 = cos(qJ(5));
t24 = sin(qJ(5));
t40 = t21 * t24;
t12 = t20 * t27 - t40;
t26 = sin(qJ(1));
t5 = t12 * t26;
t31 = t20 * t24 + t21 * t27;
t6 = t31 * t26;
t45 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t28 = cos(qJ(1));
t39 = t28 * t20;
t7 = -t27 * t39 + t28 * t40;
t8 = t31 * t28;
t44 = -t7 * r_i_i_C(1) - t8 * r_i_i_C(2);
t43 = -t31 * r_i_i_C(1) - t12 * r_i_i_C(2);
t22 = cos(qJ(2)) * pkin(2);
t42 = pkin(1) + t22 + t46;
t38 = qJ(4) * t21;
t37 = -pkin(9) - r_i_i_C(3) + pkin(8) + pkin(7);
t36 = t26 * t38 - t45;
t35 = t28 * t38 - t44;
t34 = t41 * t20;
t32 = -t43 + t46;
t30 = -sin(qJ(2)) * pkin(2) - t34;
t1 = [-t6 * r_i_i_C(1) - t5 * r_i_i_C(2) - t26 * t42 + t37 * t28, t30 * t28 + t35, -t28 * t34 + t35, t39, t44, 0; t8 * r_i_i_C(1) - t7 * r_i_i_C(2) + t37 * t26 + t28 * t42, t30 * t26 + t36, -t26 * t34 + t36, t26 * t20, t45, 0; 0, t22 + t32, t32, -t21, t43, 0;];
Ja_transl  = t1;
