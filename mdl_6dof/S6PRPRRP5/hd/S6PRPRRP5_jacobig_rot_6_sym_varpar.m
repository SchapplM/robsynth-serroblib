% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRP5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_jacobig_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:31
% EndTime: 2019-02-26 19:52:31
% DurationCPUTime: 0.09s
% Computational Cost: add. (8->8), mult. (24->19), div. (0->0), fcn. (40->8), ass. (0->14)
t93 = sin(pkin(10));
t94 = sin(pkin(6));
t105 = t93 * t94;
t95 = cos(pkin(10));
t104 = t95 * t94;
t96 = cos(pkin(6));
t98 = sin(qJ(2));
t103 = t96 * t98;
t100 = cos(qJ(2));
t102 = t93 * t100;
t101 = t95 * t100;
t99 = cos(qJ(4));
t97 = sin(qJ(4));
t1 = [0, t105, 0, -t103 * t93 + t101, t97 * t105 - (t102 * t96 + t95 * t98) * t99, 0; 0, -t104, 0, t103 * t95 + t102, -t97 * t104 - (-t101 * t96 + t93 * t98) * t99, 0; 0, t96, 0, t94 * t98, t100 * t94 * t99 + t96 * t97, 0;];
Jg_rot  = t1;
