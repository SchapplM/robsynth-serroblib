% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:12
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRPR7_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:12:10
% EndTime: 2019-02-22 10:12:10
% DurationCPUTime: 0.04s
% Computational Cost: add. (61->16), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t78 = pkin(9) + qJ(4);
t74 = sin(t78);
t79 = sin(qJ(1));
t86 = t79 * t74;
t77 = pkin(10) + qJ(6);
t75 = cos(t77);
t85 = t79 * t75;
t76 = cos(t78);
t84 = t79 * t76;
t80 = cos(qJ(1));
t83 = t80 * t74;
t82 = t80 * t75;
t81 = t80 * t76;
t73 = sin(t77);
t72 = -t79 * t73 + t74 * t82;
t71 = t73 * t83 + t85;
t70 = t80 * t73 + t74 * t85;
t69 = -t73 * t86 + t82;
t1 = [t72, 0, 0, t75 * t84, 0, t69; t70, 0, 0, -t75 * t81, 0, t71; 0, 0, 0, -t74 * t75, 0, -t76 * t73; -t71, 0, 0, -t73 * t84, 0, -t70; t69, 0, 0, t73 * t81, 0, t72; 0, 0, 0, t74 * t73, 0, -t76 * t75; -t81, 0, 0, t86, 0, 0; -t84, 0, 0, -t83, 0, 0; 0, 0, 0, t76, 0, 0;];
JR_rot  = t1;
