% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP13
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRP13_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:49
% EndTime: 2019-02-26 21:52:49
% DurationCPUTime: 0.07s
% Computational Cost: add. (8->8), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t106 = sin(pkin(6));
t110 = sin(qJ(1));
t119 = t110 * t106;
t109 = sin(qJ(2));
t118 = t110 * t109;
t112 = cos(qJ(2));
t117 = t110 * t112;
t113 = cos(qJ(1));
t116 = t113 * t106;
t115 = t113 * t109;
t114 = t113 * t112;
t111 = cos(qJ(4));
t108 = sin(qJ(4));
t107 = cos(pkin(6));
t1 = [0, t119, 0, -t107 * t118 + t114, t108 * t119 - (t107 * t117 + t115) * t111, 0; 0, -t116, 0, t107 * t115 + t117, -t108 * t116 - (-t107 * t114 + t118) * t111, 0; 1, t107, 0, t106 * t109, t106 * t111 * t112 + t107 * t108, 0;];
Jg_rot  = t1;
