% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPRP2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPRP2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPRP2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_jacobiaD_transl_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:33:01
% EndTime: 2019-02-26 19:33:02
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->24), mult. (132->34), div. (0->0), fcn. (104->4), ass. (0->19)
t82 = cos(qJ(3));
t95 = -t82 * pkin(3) - pkin(1) - pkin(2);
t94 = pkin(3) * qJD(3);
t81 = sin(qJ(1));
t93 = qJD(1) * t81;
t83 = cos(qJ(1));
t92 = qJD(1) * t83;
t91 = qJD(3) * t81;
t90 = qJD(3) * t83;
t80 = sin(qJ(3));
t89 = pkin(3) * t80 + qJ(2);
t85 = t80 * t81 + t82 * t83;
t84 = t85 * qJD(1);
t73 = -t80 * t91 - t82 * t90 + t84;
t74 = (t91 - t93) * t82 + (-t90 + t92) * t80;
t88 = -t74 * r_i_i_C(1) - t73 * r_i_i_C(2);
t87 = t73 * r_i_i_C(1) - t74 * r_i_i_C(2);
t86 = t80 * t83 - t81 * t82;
t1 = [t83 * qJD(2) + t85 * t94 + (-t89 * t81 + t95 * t83) * qJD(1) - t87, t92 (-t85 * qJD(3) + t84) * pkin(3) + t87, 0; t81 * qJD(2) - t86 * t94 + (t95 * t81 + t89 * t83) * qJD(1) - t88, t93, t88 + (-qJD(1) + qJD(3)) * pkin(3) * t86, 0; 0, 0, 0, 0;];
JaD_transl  = t1;
