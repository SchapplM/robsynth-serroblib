% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:03:30
% EndTime: 2019-02-26 20:03:31
% DurationCPUTime: 0.19s
% Computational Cost: add. (197->46), mult. (509->82), div. (0->0), fcn. (652->10), ass. (0->36)
t31 = sin(qJ(3));
t34 = cos(qJ(3));
t42 = pkin(3) + pkin(9) + r_i_i_C(2);
t54 = qJ(4) * t31 + t42 * t34 + pkin(2);
t53 = pkin(4) + pkin(8);
t52 = pkin(5) + r_i_i_C(1);
t29 = sin(pkin(6));
t51 = t29 * t31;
t50 = t29 * t34;
t30 = sin(qJ(5));
t49 = t30 * t31;
t35 = cos(qJ(2));
t48 = t30 * t35;
t33 = cos(qJ(5));
t47 = t31 * t33;
t46 = t33 * t35;
t45 = r_i_i_C(3) + qJ(6);
t44 = cos(pkin(6));
t43 = cos(pkin(10));
t28 = sin(pkin(10));
t40 = t28 * t44;
t39 = t29 * t43;
t38 = t44 * t43;
t36 = t52 * t30 - t45 * t33 + qJ(4);
t32 = sin(qJ(2));
t23 = t32 * t51 - t44 * t34;
t22 = -t32 * t40 + t43 * t35;
t21 = t43 * t32 + t35 * t40;
t20 = t28 * t35 + t32 * t38;
t19 = t28 * t32 - t35 * t38;
t15 = t23 * t33 + t29 * t48;
t13 = t22 * t31 - t28 * t50;
t11 = t20 * t31 + t34 * t39;
t3 = -t13 * t33 + t21 * t30;
t1 = -t11 * t33 + t19 * t30;
t2 = [0, t52 * (-t21 * t49 + t22 * t33) + t45 * (t21 * t47 + t22 * t30) + t53 * t22 - t54 * t21, -t42 * t13 + t36 * (t22 * t34 + t28 * t51) t13, t45 * (t13 * t30 + t21 * t33) - t52 * t3, t3; 0, t52 * (-t19 * t49 + t20 * t33) + t45 * (t19 * t47 + t20 * t30) + t53 * t20 - t54 * t19, -t42 * t11 + t36 * (t20 * t34 - t31 * t39) t11, t45 * (t11 * t30 + t19 * t33) - t52 * t1, t1; 1 (t52 * (t31 * t48 + t32 * t33) + t45 * (t30 * t32 - t31 * t46) + t53 * t32 + t54 * t35) * t29, -t42 * t23 + t36 * (t44 * t31 + t32 * t50) t23, -t45 * (-t23 * t30 + t29 * t46) + t52 * t15, -t15;];
Ja_transl  = t2;
