% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14V3_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobig_rot_6_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:09
% EndTime: 2019-04-12 15:12:09
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (24->19), div. (0->0), fcn. (42->8), ass. (0->14)
t120 = sin(qJ(2));
t121 = sin(qJ(1));
t130 = t121 * t120;
t124 = cos(qJ(2));
t129 = t121 * t124;
t119 = sin(qJ(4));
t125 = cos(qJ(1));
t128 = t125 * t119;
t127 = t125 * t120;
t123 = cos(qJ(4));
t126 = t125 * t123;
t122 = cos(qJ(5));
t118 = sin(qJ(5));
t1 = [0, t121, 0, t127, -t121 * t123 + t124 * t128 (t121 * t119 + t124 * t126) * t118 - t122 * t127; 0, -t125, 0, t130, t119 * t129 + t126 (t123 * t129 - t128) * t118 - t122 * t130; 1, 0, 0, -t124, t120 * t119, t120 * t123 * t118 + t124 * t122;];
Jg_rot  = t1;
