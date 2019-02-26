% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PPRR1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
%
% Output:
% JaD_transl [3x4]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4PPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:21:42
% EndTime: 2019-02-26 19:21:42
% DurationCPUTime: 0.02s
% Computational Cost: add. (48->10), mult. (52->17), div. (0->0), fcn. (40->6), ass. (0->14)
t68 = qJ(3) + qJ(4);
t65 = sin(t68);
t66 = cos(t68);
t67 = qJD(3) + qJD(4);
t69 = sin(pkin(6));
t70 = cos(pkin(6));
t63 = (t65 * t70 - t66 * t69) * t67;
t64 = (t65 * t69 + t66 * t70) * t67;
t75 = t63 * r_i_i_C(1) + t64 * r_i_i_C(2);
t74 = -t64 * r_i_i_C(1) + t63 * r_i_i_C(2);
t73 = pkin(3) * qJD(3);
t72 = cos(qJ(3));
t71 = sin(qJ(3));
t1 = [0, 0 (-t69 * t71 - t70 * t72) * t73 + t74, t74; 0, 0 (-t69 * t72 + t70 * t71) * t73 + t75, t75; 0, 0, 0, 0;];
JaD_transl  = t1;
