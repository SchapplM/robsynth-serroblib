% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR11_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR11_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_jacobia_transl_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:54:30
% EndTime: 2019-02-26 20:54:31
% DurationCPUTime: 0.20s
% Computational Cost: add. (228->52), mult. (532->88), div. (0->0), fcn. (694->14), ass. (0->47)
t59 = pkin(9) + sin(pkin(13)) * pkin(4);
t24 = cos(pkin(13)) * pkin(4) + pkin(3);
t27 = pkin(13) + qJ(5);
t25 = sin(t27);
t26 = cos(t27);
t58 = t26 * r_i_i_C(1) - t25 * r_i_i_C(2) + t24;
t56 = r_i_i_C(3) + pkin(10) + qJ(4);
t55 = cos(qJ(3));
t30 = sin(pkin(7));
t36 = sin(qJ(3));
t54 = t30 * t36;
t31 = sin(pkin(6));
t32 = cos(pkin(12));
t53 = t31 * t32;
t33 = cos(pkin(7));
t52 = t33 * t36;
t29 = sin(pkin(12));
t37 = sin(qJ(1));
t51 = t37 * t29;
t50 = t37 * t31;
t49 = t37 * t32;
t38 = cos(qJ(1));
t48 = t38 * t29;
t47 = t38 * t31;
t46 = t38 * t32;
t45 = t31 * qJ(2);
t44 = t30 * t55;
t43 = t33 * t55;
t42 = t31 * t44;
t34 = cos(pkin(6));
t17 = -t34 * t46 + t51;
t11 = -t17 * t30 + t33 * t47;
t40 = t34 * t49 + t48;
t39 = t40 * t33;
t18 = t34 * t48 + t49;
t4 = -t17 * t52 + t18 * t55 - t47 * t54;
t13 = t30 * t40 + t33 * t50;
t3 = t17 * t43 + t18 * t36 + t38 * t42;
t19 = -t34 * t51 + t46;
t16 = -t30 * t53 + t34 * t33;
t10 = t34 * t54 + (t55 * t29 + t32 * t52) * t31;
t9 = t31 * t29 * t36 - t34 * t44 - t43 * t53;
t8 = t19 * t55 + (t30 * t50 - t39) * t36;
t7 = t19 * t36 - t37 * t42 + t55 * t39;
t2 = t13 * t25 + t8 * t26;
t1 = t13 * t26 - t8 * t25;
t5 = [-t37 * pkin(1) - t18 * pkin(2) - t56 * t3 + t38 * t45 - t58 * t4 + (t25 * r_i_i_C(1) + t26 * r_i_i_C(2) + t59) * t11, t50, t56 * t8 - t58 * t7, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2), 0; t38 * pkin(1) + t19 * pkin(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t59 * t13 + t8 * t24 + t37 * t45 + t56 * t7, -t47, -t3 * t58 + t56 * t4, t3 (-t11 * t26 - t4 * t25) * r_i_i_C(1) + (t11 * t25 - t4 * t26) * r_i_i_C(2), 0; 0, t34, t56 * t10 - t58 * t9, t9 (-t10 * t25 + t16 * t26) * r_i_i_C(1) + (-t10 * t26 - t16 * t25) * r_i_i_C(2), 0;];
Ja_transl  = t5;
