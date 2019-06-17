% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-06-12 14:37
% Revision: aab8d7cd0cba739f5e0ec8d53b8419901d1154b0 (2019-06-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RPRRR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobigD_rot_5_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-06-12 14:37:34
% EndTime: 2019-06-12 14:37:34
% DurationCPUTime: 0.05s
% Computational Cost: add. (11->9), mult. (41->23), div. (0->0), fcn. (41->6), ass. (0->14)
t72 = sin(qJ(1));
t82 = qJD(1) * t72;
t75 = cos(qJ(1));
t81 = qJD(1) * t75;
t71 = sin(qJ(3));
t80 = qJD(3) * t71;
t74 = cos(qJ(3));
t79 = qJD(3) * t74;
t78 = qJD(3) * t75;
t77 = qJD(1) * t74 - qJD(4);
t73 = cos(qJ(4));
t76 = (qJD(4) * t74 - qJD(1)) * t73;
t70 = sin(qJ(4));
t1 = [0, 0, t81, -t71 * t82 + t74 * t78, t75 * t76 + (-t71 * t78 - t77 * t72) * t70; 0, 0, t82, t71 * t81 + t72 * t79, t72 * t76 + (-t72 * t80 + t77 * t75) * t70; 0, 0, 0, t80, t71 * qJD(4) * t73 + t70 * t79;];
JgD_rot  = t1;
