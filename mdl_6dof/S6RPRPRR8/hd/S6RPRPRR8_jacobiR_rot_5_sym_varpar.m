% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRR8_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:34:42
% EndTime: 2019-02-22 10:34:42
% DurationCPUTime: 0.05s
% Computational Cost: add. (37->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t74 = sin(qJ(5));
t75 = sin(qJ(1));
t81 = t75 * t74;
t76 = cos(qJ(5));
t80 = t75 * t76;
t77 = cos(qJ(1));
t79 = t77 * t74;
t78 = t77 * t76;
t73 = qJ(3) + pkin(10);
t72 = cos(t73);
t71 = sin(t73);
t70 = t71 * t78 - t81;
t69 = t71 * t79 + t80;
t68 = t71 * t80 + t79;
t67 = -t71 * t81 + t78;
t1 = [t70, 0, t72 * t80, 0, t67, 0; t68, 0, -t72 * t78, 0, t69, 0; 0, 0, -t71 * t76, 0, -t72 * t74, 0; -t69, 0, -t72 * t81, 0, -t68, 0; t67, 0, t72 * t79, 0, t70, 0; 0, 0, t71 * t74, 0, -t72 * t76, 0; -t77 * t72, 0, t75 * t71, 0, 0, 0; -t75 * t72, 0, -t77 * t71, 0, 0, 0; 0, 0, t72, 0, 0, 0;];
JR_rot  = t1;
