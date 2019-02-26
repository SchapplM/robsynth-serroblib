% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP10_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP10_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:14:19
% EndTime: 2019-02-26 22:14:19
% DurationCPUTime: 0.18s
% Computational Cost: add. (218->47), mult. (427->80), div. (0->0), fcn. (540->12), ass. (0->36)
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t19 = cos(pkin(11)) * pkin(4) + pkin(3);
t22 = pkin(11) + qJ(5);
t20 = sin(t22);
t21 = cos(t22);
t33 = t21 * r_i_i_C(1) - t20 * r_i_i_C(2) + t19;
t43 = r_i_i_C(3) + pkin(10) + qJ(4);
t44 = t43 * t26 + t33 * t29 + pkin(2);
t42 = cos(qJ(1));
t24 = sin(pkin(6));
t28 = sin(qJ(1));
t41 = t24 * t28;
t40 = t24 * t29;
t30 = cos(qJ(2));
t39 = t24 * t30;
t38 = cos(pkin(6));
t37 = sin(pkin(11)) * pkin(4) + pkin(9);
t36 = t24 * t42;
t27 = sin(qJ(2));
t34 = t38 * t42;
t12 = t27 * t34 + t28 * t30;
t4 = t12 * t29 - t26 * t36;
t35 = t28 * t38;
t3 = t12 * t26 + t29 * t36;
t32 = t20 * r_i_i_C(1) + t21 * r_i_i_C(2) + t37;
t14 = -t27 * t35 + t42 * t30;
t13 = t42 * t27 + t30 * t35;
t11 = t28 * t27 - t30 * t34;
t10 = t38 * t26 + t27 * t40;
t9 = t24 * t27 * t26 - t38 * t29;
t8 = t14 * t29 + t26 * t41;
t7 = t14 * t26 - t28 * t40;
t2 = t13 * t20 + t8 * t21;
t1 = t13 * t21 - t8 * t20;
t5 = [-t28 * pkin(1) - t12 * pkin(2) + pkin(8) * t36 - t32 * t11 - t43 * t3 - t33 * t4, -t13 * t44 + t32 * t14, -t33 * t7 + t43 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t42 * pkin(1) + t14 * pkin(2) + pkin(8) * t41 + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t37 * t13 + t8 * t19 + t43 * t7, -t11 * t44 + t32 * t12, -t33 * t3 + t43 * t4, t3 (t11 * t21 - t4 * t20) * r_i_i_C(1) + (-t11 * t20 - t4 * t21) * r_i_i_C(2), 0; 0 (t32 * t27 + t44 * t30) * t24, t43 * t10 - t33 * t9, t9 (-t10 * t20 - t21 * t39) * r_i_i_C(1) + (-t10 * t21 + t20 * t39) * r_i_i_C(2), 0;];
Ja_transl  = t5;
