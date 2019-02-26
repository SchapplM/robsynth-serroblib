% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:04:21
% EndTime: 2019-02-26 21:04:21
% DurationCPUTime: 0.12s
% Computational Cost: add. (154->35), mult. (158->46), div. (0->0), fcn. (103->8), ass. (0->32)
t56 = qJ(3) + qJ(4);
t50 = pkin(10) + t56;
t48 = sin(t50);
t49 = cos(t50);
t82 = r_i_i_C(1) * t48 + r_i_i_C(2) * t49;
t53 = cos(t56);
t74 = r_i_i_C(2) * t48;
t80 = pkin(4) * t53 - t74;
t52 = sin(t56);
t81 = -pkin(4) * t52 - t82;
t55 = qJD(3) + qJD(4);
t77 = pkin(4) * t55;
t75 = r_i_i_C(1) * t49;
t60 = cos(qJ(1));
t72 = t55 * t60;
t71 = pkin(3) * qJD(3);
t58 = sin(qJ(1));
t70 = qJD(1) * t58;
t51 = qJD(1) * t60;
t69 = -pkin(1) - r_i_i_C(3) - qJ(5) - pkin(8) - pkin(7);
t67 = qJD(1) * t75;
t68 = t58 * t67 + t82 * t72;
t59 = cos(qJ(3));
t66 = t59 * t71;
t64 = qJD(1) * (pkin(3) * t59 + t80);
t57 = sin(qJ(3));
t63 = pkin(3) * t57 + qJ(2) - t81;
t62 = (-t75 - t80) * t55;
t61 = qJD(2) + t53 * t77 + t66 + (-t74 + t75) * t55;
t44 = t60 * t67;
t39 = -t52 * t77 - t57 * t71;
t1 = [-t58 * qJD(5) + t61 * t60 + (-t63 * t58 + t69 * t60) * qJD(1), t51, t44 + t60 * t64 + (-t55 * t82 + t39) * t58, t81 * t58 * t55 + t80 * t51 + t44, -t70, 0; qJD(5) * t60 + t61 * t58 + (t69 * t58 + t63 * t60) * qJD(1), t70, -t60 * t39 + t58 * t64 + t68, -t70 * t74 + (t52 * t72 + t53 * t70) * pkin(4) + t68, t51, 0; 0, 0, t62 - t66, t62, 0, 0;];
JaD_transl  = t1;
