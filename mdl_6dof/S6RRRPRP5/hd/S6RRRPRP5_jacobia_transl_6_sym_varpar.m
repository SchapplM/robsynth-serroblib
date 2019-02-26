% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:11:36
% EndTime: 2019-02-26 22:11:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (271->35), mult. (217->46), div. (0->0), fcn. (249->10), ass. (0->31)
t43 = pkin(5) + r_i_i_C(1);
t36 = r_i_i_C(3) + qJ(6);
t24 = qJ(3) + pkin(10);
t13 = pkin(4) * cos(t24) + cos(qJ(3)) * pkin(3);
t11 = pkin(2) + t13;
t27 = cos(qJ(2));
t25 = sin(qJ(2));
t41 = r_i_i_C(2) + pkin(9) + qJ(4) + pkin(8);
t33 = t41 * t25;
t45 = t27 * t11 + pkin(1) + t33;
t21 = qJ(5) + t24;
t19 = sin(t21);
t20 = cos(t21);
t44 = t36 * t19 + t43 * t20 + t11;
t12 = pkin(4) * sin(t24) + sin(qJ(3)) * pkin(3);
t42 = pkin(7) + t12;
t26 = sin(qJ(1));
t39 = t26 * t27;
t28 = cos(qJ(1));
t38 = t28 * t19;
t37 = t28 * t20;
t35 = t36 * t20 * t25;
t34 = t43 * t19;
t7 = t19 * t39 + t37;
t8 = t20 * t39 - t38;
t31 = t36 * t8 - t43 * t7;
t10 = t26 * t19 + t27 * t37;
t9 = -t26 * t20 + t27 * t38;
t30 = t36 * t10 - t43 * t9;
t29 = -t44 * t25 + t41 * t27;
t1 = [-t45 * t26 + t42 * t28 - t36 * t7 - t43 * t8, t29 * t28, -t28 * t27 * t12 + t26 * t13 + t30, t28 * t25, t30, t9; t43 * t10 + t42 * t26 + t45 * t28 + t36 * t9, t29 * t26, -t12 * t39 - t28 * t13 + t31, t26 * t25, t31, t7; 0, t44 * t27 + t33 (-t12 - t34) * t25 + t35, -t27, -t25 * t34 + t35, t25 * t19;];
Ja_transl  = t1;
