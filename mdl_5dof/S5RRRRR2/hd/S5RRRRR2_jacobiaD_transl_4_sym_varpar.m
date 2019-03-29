% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_jacobiaD_transl_4_sym_varpar: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:51
% EndTime: 2019-03-29 15:26:51
% DurationCPUTime: 0.14s
% Computational Cost: add. (175->28), mult. (154->44), div. (0->0), fcn. (99->8), ass. (0->33)
t57 = qJ(3) + qJ(4);
t51 = sin(t57);
t58 = qJ(1) + qJ(2);
t52 = sin(t58);
t54 = cos(t58);
t56 = qJD(1) + qJD(2);
t77 = t54 * t56;
t53 = cos(t57);
t55 = qJD(3) + qJD(4);
t78 = t53 * t55;
t83 = t51 * t77 + t52 * t78;
t59 = sin(qJ(3));
t82 = pkin(2) * t59;
t81 = r_i_i_C(2) * t53;
t80 = t51 * t55;
t79 = t52 * t56;
t76 = pkin(1) * qJD(1);
t60 = cos(qJ(3));
t75 = qJD(3) * t60;
t74 = r_i_i_C(1) * t78;
t73 = t56 * t81;
t72 = t52 * t80;
t71 = t51 * t79;
t68 = qJD(3) * t82;
t67 = -pkin(2) * t60 - r_i_i_C(1) * t53;
t66 = -r_i_i_C(1) * t51 - t81;
t65 = t66 * t55;
t64 = r_i_i_C(1) * t71 + t52 * t73 + (r_i_i_C(2) * t80 - t74) * t54;
t63 = t65 - t68;
t62 = r_i_i_C(1) * t72 + t52 * t68 + (-r_i_i_C(3) * t52 + t67 * t54) * t56 + t83 * r_i_i_C(2);
t61 = r_i_i_C(2) * t71 + r_i_i_C(3) * t77 + t63 * t54 + t67 * t79;
t41 = r_i_i_C(2) * t72;
t1 = [-cos(qJ(1)) * t76 + t62, t62 (-t54 * t75 + t59 * t79) * pkin(2) + t64, t64, 0; -sin(qJ(1)) * t76 + t61, t61, t41 + (-pkin(2) * t75 - t74) * t52 + (t66 - t82) * t77, -t83 * r_i_i_C(1) - t54 * t73 + t41, 0; 0, 0, t63, t65, 0;];
JaD_transl  = t1;
