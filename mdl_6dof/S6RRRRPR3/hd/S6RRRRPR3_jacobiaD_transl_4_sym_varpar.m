% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:52
% EndTime: 2019-02-26 22:31:52
% DurationCPUTime: 0.14s
% Computational Cost: add. (202->36), mult. (170->47), div. (0->0), fcn. (109->8), ass. (0->37)
t54 = sin(qJ(2));
t51 = qJD(2) + qJD(3);
t47 = qJD(4) + t51;
t53 = qJ(2) + qJ(3);
t50 = qJ(4) + t53;
t46 = cos(t50);
t77 = r_i_i_C(2) * t46;
t45 = sin(t50);
t79 = r_i_i_C(1) * t45;
t62 = t77 + t79;
t60 = t62 * t47;
t48 = sin(t53);
t80 = pkin(3) * t48;
t58 = -t51 * t80 - t60;
t73 = pkin(2) * qJD(2);
t82 = -t54 * t73 + t58;
t75 = t46 * t47;
t68 = r_i_i_C(1) * t75;
t49 = cos(t53);
t74 = t49 * t51;
t81 = -pkin(3) * t74 - t68;
t78 = r_i_i_C(2) * t45;
t76 = r_i_i_C(3) + pkin(9) + pkin(8) + pkin(7);
t55 = sin(qJ(1));
t72 = qJD(1) * t55;
t57 = cos(qJ(1));
t71 = qJD(1) * t57;
t67 = t47 * t78;
t65 = qJD(1) * t77;
t66 = t55 * t65 + t57 * t67 + t72 * t79;
t56 = cos(qJ(2));
t63 = -t56 * t73 + t81;
t61 = -pkin(2) * t56 - pkin(3) * t49 - r_i_i_C(1) * t46 - pkin(1) + t78;
t59 = -t57 * t68 + t66;
t44 = -pkin(2) * t54 - t80;
t39 = t55 * t67;
t1 = [-t82 * t55 + (-t76 * t55 + t61 * t57) * qJD(1), -t44 * t72 + t63 * t57 + t66 (t48 * t72 - t57 * t74) * pkin(3) + t59, t59, 0, 0; t82 * t57 + (t61 * t55 + t76 * t57) * qJD(1), t39 + t63 * t55 + (t44 - t62) * t71, t39 + t81 * t55 + (-t62 - t80) * t71, -t57 * t65 + t39 + (-t45 * t71 - t55 * t75) * r_i_i_C(1), 0, 0; 0, t82, t58, -t60, 0, 0;];
JaD_transl  = t1;
