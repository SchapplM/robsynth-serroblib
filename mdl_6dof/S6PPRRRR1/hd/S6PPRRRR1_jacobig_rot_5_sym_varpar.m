% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR1_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobig_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:43
% DurationCPUTime: 0.04s
% Computational Cost: add. (25->13), mult. (77->32), div. (0->0), fcn. (111->10), ass. (0->20)
t108 = sin(pkin(12));
t114 = cos(pkin(6));
t120 = t108 * t114;
t109 = sin(pkin(7));
t110 = sin(pkin(6));
t119 = t109 * t110;
t113 = cos(pkin(7));
t118 = t110 * t113;
t112 = cos(pkin(12));
t117 = t112 * t114;
t116 = cos(qJ(3));
t115 = sin(qJ(3));
t111 = cos(pkin(13));
t107 = sin(pkin(13));
t106 = -t112 * t107 - t111 * t120;
t105 = -t108 * t107 + t111 * t117;
t104 = -t114 * t109 * t116 + (-t111 * t113 * t116 + t107 * t115) * t110;
t103 = (-t107 * t120 + t112 * t111) * t115 + (-t106 * t113 - t108 * t119) * t116;
t102 = (t107 * t117 + t108 * t111) * t115 + (-t105 * t113 + t112 * t119) * t116;
t1 = [0, 0, -t106 * t109 + t108 * t118, t103, t103, 0; 0, 0, -t105 * t109 - t112 * t118, t102, t102, 0; 0, 0, -t111 * t119 + t114 * t113, t104, t104, 0;];
Jg_rot  = t1;
