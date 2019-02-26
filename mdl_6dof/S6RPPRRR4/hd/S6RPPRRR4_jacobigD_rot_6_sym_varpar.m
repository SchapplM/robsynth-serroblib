% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RPPRRR4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobigD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.02s
% Computational Cost: add. (16->7), mult. (46->12), div. (0->0), fcn. (50->6), ass. (0->14)
t128 = sin(qJ(4));
t135 = qJD(4) * t128;
t134 = qJD(4) * cos(qJ(4));
t126 = sin(pkin(10));
t127 = cos(pkin(10));
t129 = sin(qJ(1));
t131 = cos(qJ(1));
t133 = t131 * t126 - t129 * t127;
t132 = t129 * t126 + t131 * t127;
t125 = t133 * qJD(1);
t124 = t132 * qJD(1);
t123 = t125 * t128 + t132 * t134;
t122 = t124 * t128 - t133 * t134;
t1 = [0, 0, 0, -t124, t123, t123; 0, 0, 0, t125, t122, t122; 0, 0, 0, 0, -t135, -t135;];
JgD_rot  = t1;
