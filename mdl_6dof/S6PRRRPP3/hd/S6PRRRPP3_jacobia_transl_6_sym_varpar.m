% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPP3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:01
% EndTime: 2019-02-26 20:10:01
% DurationCPUTime: 0.17s
% Computational Cost: add. (228->44), mult. (590->82), div. (0->0), fcn. (758->10), ass. (0->38)
t29 = sin(qJ(3));
t32 = cos(qJ(3));
t41 = pkin(5) + pkin(9) + r_i_i_C(1);
t51 = pkin(3) * t32 + t41 * t29 + pkin(2);
t27 = sin(pkin(6));
t50 = t27 * t29;
t49 = t27 * t32;
t33 = cos(qJ(2));
t48 = t27 * t33;
t28 = sin(qJ(4));
t47 = t28 * t32;
t31 = cos(qJ(4));
t46 = t31 * t32;
t45 = t32 * t33;
t44 = r_i_i_C(2) + qJ(5);
t43 = cos(pkin(6));
t42 = cos(pkin(10));
t40 = pkin(4) + r_i_i_C(3) + qJ(6);
t26 = sin(pkin(10));
t38 = t26 * t43;
t37 = t27 * t42;
t36 = t43 * t42;
t34 = t44 * t28 + t40 * t31 + pkin(3);
t30 = sin(qJ(2));
t24 = t43 * t29 + t30 * t49;
t22 = -t30 * t38 + t42 * t33;
t21 = t42 * t30 + t33 * t38;
t20 = t26 * t33 + t30 * t36;
t19 = t26 * t30 - t33 * t36;
t14 = t24 * t31 - t28 * t48;
t13 = t24 * t28 + t31 * t48;
t12 = t22 * t32 + t26 * t50;
t10 = t20 * t32 - t29 * t37;
t4 = t12 * t31 + t21 * t28;
t3 = t12 * t28 - t21 * t31;
t2 = t10 * t31 + t19 * t28;
t1 = t10 * t28 - t19 * t31;
t5 = [0, t22 * pkin(8) + t44 * (-t21 * t47 - t22 * t31) + t40 * (-t21 * t46 + t22 * t28) - t51 * t21, t41 * t12 + t34 * (-t22 * t29 + t26 * t49) -t40 * t3 + t44 * t4, t3, t4; 0, t20 * pkin(8) + t44 * (-t19 * t47 - t20 * t31) + t40 * (-t19 * t46 + t20 * t28) - t51 * t19, t41 * t10 + t34 * (-t20 * t29 - t32 * t37) -t40 * t1 + t44 * t2, t1, t2; 1 (t44 * (t28 * t45 - t30 * t31) + t40 * (t28 * t30 + t31 * t45) + pkin(8) * t30 + t51 * t33) * t27, t41 * t24 + t34 * (-t30 * t50 + t43 * t32) -t40 * t13 + t44 * t14, t13, t14;];
Ja_transl  = t5;
