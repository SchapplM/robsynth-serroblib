% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:51:15
% EndTime: 2019-02-26 20:51:15
% DurationCPUTime: 0.21s
% Computational Cost: add. (253->43), mult. (382->70), div. (0->0), fcn. (331->7), ass. (0->34)
t67 = pkin(10) + qJ(3);
t65 = sin(t67);
t66 = cos(t67);
t69 = sin(qJ(5));
t71 = cos(qJ(5));
t103 = -t65 * t71 + t66 * t69;
t101 = qJD(3) - qJD(5);
t78 = t65 * t69 + t66 * t71;
t59 = t101 * t78;
t98 = pkin(3) + pkin(4);
t102 = (-qJ(4) * t66 + t98 * t65) * qJD(3) - t65 * qJD(4);
t92 = qJD(3) * t65;
t100 = t103 * qJD(5) + t71 * t92;
t70 = sin(qJ(1));
t94 = qJD(1) * t70;
t72 = cos(qJ(1));
t93 = qJD(1) * t72;
t91 = qJD(3) * t66;
t89 = qJ(4) * qJD(1);
t88 = qJ(4) * qJD(3);
t87 = pkin(8) + r_i_i_C(3) - pkin(7) - qJ(2);
t83 = t72 * t91;
t76 = qJD(1) * t103;
t54 = t59 * t72 + t70 * t76;
t75 = qJD(1) * t78;
t55 = t100 * t72 - t69 * t83 + t70 * t75;
t82 = t54 * r_i_i_C(1) + t55 * r_i_i_C(2);
t56 = -t70 * t59 + t72 * t76;
t57 = t72 * t75 + (t69 * t91 - t100) * t70;
t81 = -t56 * r_i_i_C(1) - t57 * r_i_i_C(2);
t80 = -r_i_i_C(1) * t101 * t103 - t59 * r_i_i_C(2);
t77 = -t98 * qJD(3) + qJD(4);
t74 = -qJ(4) * t65 - t98 * t66 - cos(pkin(10)) * pkin(2) - pkin(1);
t1 = [-t57 * r_i_i_C(1) + t56 * r_i_i_C(2) + t72 * qJD(2) + t102 * t70 + (t87 * t70 + t74 * t72) * qJD(1), t93 (-t72 * t88 + t98 * t94) * t65 + (-t70 * t89 + t77 * t72) * t66 - t82, -t65 * t94 + t83, t82, 0; -t55 * r_i_i_C(1) + t54 * r_i_i_C(2) + t70 * qJD(2) - t102 * t72 + (t74 * t70 - t87 * t72) * qJD(1), t94 (-t70 * t88 - t98 * t93) * t65 + (t77 * t70 + t72 * t89) * t66 - t81, t65 * t93 + t70 * t91, t81, 0; 0, 0, -t102 - t80, t92, t80, 0;];
JaD_transl  = t1;
