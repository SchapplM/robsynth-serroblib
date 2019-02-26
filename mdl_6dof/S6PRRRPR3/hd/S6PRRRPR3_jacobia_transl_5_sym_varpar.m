% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:42
% EndTime: 2019-02-26 20:11:42
% DurationCPUTime: 0.15s
% Computational Cost: add. (175->29), mult. (265->51), div. (0->0), fcn. (329->10), ass. (0->30)
t46 = pkin(4) - r_i_i_C(2);
t45 = r_i_i_C(3) + qJ(5);
t44 = r_i_i_C(1) + pkin(9) + pkin(8);
t25 = sin(pkin(11));
t26 = sin(pkin(6));
t43 = t25 * t26;
t27 = cos(pkin(11));
t42 = t26 * t27;
t30 = sin(qJ(2));
t41 = t26 * t30;
t31 = cos(qJ(3));
t40 = t26 * t31;
t28 = cos(pkin(6));
t39 = t28 * t30;
t32 = cos(qJ(2));
t38 = t28 * t32;
t17 = t25 * t32 + t27 * t39;
t24 = qJ(3) + qJ(4);
t22 = sin(t24);
t23 = cos(t24);
t7 = t17 * t22 + t23 * t42;
t37 = t45 * (t17 * t23 - t22 * t42) - t46 * t7;
t19 = -t25 * t39 + t27 * t32;
t9 = t19 * t22 - t23 * t43;
t36 = -t46 * t9 + t45 * (t19 * t23 + t22 * t43);
t14 = t22 * t41 - t28 * t23;
t35 = t45 * (t28 * t22 + t23 * t41) - t46 * t14;
t34 = t31 * pkin(3) + t45 * t22 + t46 * t23 + pkin(2);
t29 = sin(qJ(3));
t1 = [0, t44 * t19 + t34 * (-t25 * t38 - t27 * t30) (-t19 * t29 + t25 * t40) * pkin(3) + t36, t36, t9, 0; 0, t44 * t17 + t34 * (-t25 * t30 + t27 * t38) (-t17 * t29 - t27 * t40) * pkin(3) + t37, t37, t7, 0; 1 (t44 * t30 + t34 * t32) * t26 (t28 * t31 - t29 * t41) * pkin(3) + t35, t35, t14, 0;];
Ja_transl  = t1;
